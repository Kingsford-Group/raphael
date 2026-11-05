namespace QuantumSketchExperiments {
    open Std.Canon;
    open Std.Measurement;
    open Std.Random;

    // Small-instance functional validation of coherent score-based search.
    // The search never receives a precomputed list of scores. Each marking
    // circuit computes QAE coherently from the underlying input sets.
    // The tractable simulator preset uses ONE coherent QAE register, rather
    // than the logarithmic coherent-median amplification used in the theorem.
    // Measured classical thresholds use the median of three separate QAE runs.

    // On the QPE basis state |v>, p=sin^2(pi*v/2^phaseBits).
    // This classical enumeration synthesizes the reversible numerical decoder
    // and comparison. It does NOT enumerate input intersections or scores.
    operation SearchMarkScore(db : Int[][], target : Int[], pairMode : Bool,
                              threshold : Double, margin : Double,
                              candidate : Qubit[], phase : Qubit[], good : Qubit) : Unit is Adj + Ctl {
        let bits = BitsNeeded(Length(db));
        let phaseBits = Length(phase);
        let resolution = 1 <<< phaseBits;
        if pairMode {
            let candidateCapacity = 1 <<< (2 * bits);
            for i in 0..Length(db) - 1 {
                for j in i + 1..Length(db) - 1 {
                    let encoded = i + (j <<< bits);
                    for value in 0..resolution - 1 {
                        let p = PhaseEstimateProbability(value, phaseBits);
                        let score = ContainmentToJaccard(p, Length(db[i]), Length(db[j]));
                        if score > threshold + margin {
                            ApplyControlledOnInt(encoded + value * candidateCapacity,
                                                 X, candidate + phase, good);
                        }
                    }
                }
            }
        } else {
            for value in 0..resolution - 1 {
                if PhaseEstimateProbability(value, phaseBits) > threshold + margin {
                    ApplyControlledOnInt(value, X, phase, good);
                }
            }
        }
    }

    operation SearchPrepareMarked(db : Int[][], target : Int[], pairMode : Bool,
                                  threshold : Double, margin : Double,
                                  candidate : Qubit[], sample : Qubit[], key : Qubit[],
                                  membership : Qubit, phase : Qubit[], good : Qubit) : Unit is Adj + Ctl {
        let bits = BitsNeeded(Length(db));
        if pairMode {
            PrepareIndex(Length(db), candidate[0..bits - 1]);
            PrepareIndex(Length(db), candidate[bits..2 * bits - 1]);
        } else {
            PrepareIndex(Length(db), candidate);
        }
        CoherentCandidateQae(db, target, pairMode, candidate, sample, key, membership, phase);
        SearchMarkScore(db, target, pairMode, threshold, margin, candidate, phase, good);
    }

    // Return (candidate, success flag, logical A/A-dagger calls).
    // The inverse and all QAE work registers participate in Grover diffusion.
    operation SearchAmplifiedCandidate(db : Int[][], target : Int[], pairMode : Bool,
                                       phaseBits : Int, threshold : Double,
                                       margin : Double, iterations : Int) : (Int, Bool, Int) {
        let bits = BitsNeeded(Length(db));
        let candidateBits = pairMode ? 2 * bits | bits;
        use candidate = Qubit[candidateBits];
        use sample = Qubit[DataSampleWidth(db)];
        use key = Qubit[DataKeyWidth(db, target)];
        use membership = Qubit();
        use phase = Qubit[phaseBits];
        use good = Qubit();
        SearchPrepareMarked(db, target, pairMode, threshold, margin,
                            candidate, sample, key, membership, phase, good);
        for _ in 1..iterations {
            Z(good);
            Adjoint SearchPrepareMarked(db, target, pairMode, threshold, margin,
                                        candidate, sample, key, membership, phase, good);
            ReflectZero(candidate + sample + key + [membership] + phase + [good]);
            SearchPrepareMarked(db, target, pairMode, threshold, margin,
                                candidate, sample, key, membership, phase, good);
        }
        let succeeded = MResetZ(good) == One;
        let encoded = MeasureInteger(candidate);
        ResetAll(sample + key + [membership] + phase);
        return (encoded, succeeded, (1 + 2 * iterations) * (2 * (1 <<< phaseBits) - 1));
    }

    operation SearchThreshold(db : Int[][], target : Int[], pairMode : Bool,
                              phaseBits : Int, encoded : Int) : (Double, Int) {
        return CandidateEstimate(db, target, pairMode, encoded, phaseBits, 3);
    }

    // Return (encoded candidate, measured threshold, logical A/A-dagger calls).
    // For pairs, encoded=i+2^BitsNeeded(R)*j, and only i<j may be returned.
    // margin is the marking margin itself, not the theorem's epsilon.
    // roundsPerUpdate and maxUpdates are explicit finite-size search budgets.
    operation SearchDatabase(db : Int[][], target : Int[], pairMode : Bool,
                             phaseBits : Int, margin : Double,
                             roundsPerUpdate : Int, maxUpdates : Int) : (Int, Double, Int) {
        if Length(db) == 0 { fail "The database must be nonempty."; }
        if pairMode and Length(db) < 2 { fail "Best-pair search requires R >= 2."; }
        if phaseBits < 1 or roundsPerUpdate < 1 or maxUpdates < 1 or margin < 0.0 {
            fail "Invalid search precision, margin, or budget.";
        }
        for source in db { ValidateSet(source, false); }
        ValidateSet(target, true);
        let bits = BitsNeeded(Length(db));
        mutable current = pairMode ? (1 <<< bits) | 0;
        let (initialThreshold, initialCalls) = SearchThreshold(db, target, pairMode, phaseBits, current);
        mutable threshold = initialThreshold;
        mutable calls = initialCalls;
        let count = pairMode ? Length(db) * Length(db) | Length(db);
        mutable cap = 1;
        while cap * cap < count { set cap += 1; }
        mutable running = true;
        mutable update = 0;
        while running and update < maxUpdates {
            mutable found = false;
            mutable round = 0;
            while not found and round < roundsPerUpdate {
                let iterations = DrawRandomInt(0, cap - 1);
                let (proposed, flagged, used) = SearchAmplifiedCandidate(
                    db, target, pairMode, phaseBits, threshold, margin, iterations);
                set calls += used;
                if flagged {
                    let (verified, verificationCalls) = SearchThreshold(
                        db, target, pairMode, phaseBits, proposed);
                    set calls += verificationCalls;
                    // A finite-size verification guard reduces false updates
                    // from the unamplified coherent estimator. Its own failures
                    // remain visible in the trial-level output and success rate.
                    if verified > threshold + margin / 2.0 {
                        set current = proposed;
                        set threshold = verified;
                        set found = true;
                    }
                }
                set round += 1;
            }
            set update += 1;
            if not found { set running = false; }
        }
        return (current, threshold, calls);
    }
}
