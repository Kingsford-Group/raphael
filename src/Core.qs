namespace QuantumSketchExperiments {
    open Std.Canon;
    open Std.Math;
    open Std.Convert;
    open Std.Measurement;
    open Std.StatePreparation;

    function BitsNeeded(n : Int) : Int {
        mutable width = 1;
        while (1 <<< width) < n { set width += 1; }
        return width;
    }

    function ValidateSet(keys : Int[], allowEmpty : Bool) : Unit {
        if not allowEmpty and Length(keys) == 0 { fail "Source set must be nonempty."; }
        for i in 0..Length(keys)-1 {
            if keys[i] < 0 { fail "Keys must be nonnegative."; }
            for j in 0..i-1 {
                if keys[i] == keys[j] { fail "Pass distinct set elements, not repeated k-mer occurrences."; }
            }
        }
    }

    function DataKeyWidth(db : Int[][], target : Int[]) : Int {
        mutable maximum = 0;
        for row in db { for x in row { if x > maximum { set maximum = x; } } }
        for x in target { if x > maximum { set maximum = x; } }
        return BitsNeeded(maximum + 1);
    }

    function DataSampleWidth(db : Int[][]) : Int {
        mutable maximum = 1;
        for row in db { if Length(row) > maximum { set maximum = Length(row); } }
        return BitsNeeded(maximum);
    }

    operation XorConstant(value : Int, register : Qubit[]) : Unit is Adj + Ctl {
        for bit in 0..Length(register)-1 {
            if ((value >>> bit) &&& 1) == 1 { X(register[bit]); }
        }
    }

    // Ancilla-free full unitary, including on inputs outside its preparation
    // subspace. This matters inside A^dagger, controlled Q, and amplitude
    // amplification: a preparation routine with input-dependent ancilla cleanup
    // is insufficient there for non-power-of-two cardinalities.
    operation PrepareIndex(nStates : Int, register : Qubit[]) : Unit is Adj + Ctl {
        let width = Length(register);
        for bit in width-1..-1..0 {
            let blockSize = 1 <<< (bit + 1);
            let half = 1 <<< bit;
            for prefix in 0..(1 <<< (width - bit - 1))-1 {
                let remaining = nStates - prefix * blockSize;
                let total = MinI(blockSize, MaxI(0, remaining));
                if total > 0 {
                    let left = MinI(half, total);
                    let angle = 2.0 * ArcCos(Sqrt(IntAsDouble(left) / IntAsDouble(total)));
                    if bit == width - 1 { Ry(angle, register[bit]); }
                    else {
                        ApplyControlledOnInt(prefix, Ry, register[bit+1..width-1], (angle, register[bit]));
                    }
                }
            }
        }
    }

    // Explicit reversible multiplexing: a small-circuit emulator, NOT physical QRAM.
    // Only input keys are embedded; intersections and similarity scores are not.
    operation XorLookup(data : Int[], address : Qubit[], output : Qubit[]) : Unit is Adj + Ctl {
        for i in 0..Length(data)-1 {
            ApplyControlledOnInt(i, XorConstant, address, (data[i], output));
        }
    }

    operation Membership(target : Int[], key : Qubit[], flag : Qubit) : Unit is Adj + Ctl {
        for x in target { ApplyControlledOnInt(x, X, key, flag); }
    }

    operation ReflectZero(register : Qubit[]) : Unit is Adj + Ctl {
        within { for q in register { X(q); } }
        apply {
            if Length(register) == 1 { Z(register[0]); }
            else { Controlled Z(register[0..Length(register)-2], register[Length(register)-1]); }
        }
    }

    operation GlobalMinus() : Unit is Adj + Ctl {
        use q = Qubit();
        X(q); Z(q); X(q);
    }

    // Candidate layout: low rbits = source i; next rbits = target j for pair search.
    // Retaining the sampled index avoids needing a reverse key-to-index dictionary.
    operation PrepareCandidateExperiment(
        db : Int[][], target : Int[], pairMode : Bool, candidate : Qubit[],
        sample : Qubit[], key : Qubit[], flag : Qubit
    ) : Unit is Adj + Ctl {
        let rbits = BitsNeeded(Length(db));
        let source = candidate[0..rbits-1];
        for i in 0..Length(db)-1 {
            ApplyControlledOnInt(i, PrepareIndex, source, (Length(db[i]), sample));
            ApplyControlledOnInt(i, XorLookup, source, (db[i], sample, key));
        }
        if pairMode {
            let targetIndex = candidate[rbits..2*rbits-1];
            for j in 0..Length(db)-1 {
                ApplyControlledOnInt(j, Membership, targetIndex, (db[j], key, flag));
            }
        } else { Membership(target, key, flag); }
    }

    // Q = - A S_0 A^dagger S_good. The minus sign matters under phase estimation.
    operation CandidateGrover(
        db : Int[][], target : Int[], pairMode : Bool, candidate : Qubit[],
        sample : Qubit[], key : Qubit[], flag : Qubit
    ) : Unit is Adj + Ctl {
        Z(flag);
        Adjoint PrepareCandidateExperiment(db, target, pairMode, candidate, sample, key, flag);
        ReflectZero(sample + key + [flag]);
        PrepareCandidateExperiment(db, target, pairMode, candidate, sample, key, flag);
        GlobalMinus();
    }

    operation FullQft(register : Qubit[]) : Unit is Adj + Ctl {
        ApplyQFT(register);
        SwapReverseRegister(register);
    }

    // Canonical, coherent phase-estimation QAE. No measurements occur here.
    operation CoherentCandidateQae(
        db : Int[][], target : Int[], pairMode : Bool, candidate : Qubit[],
        sample : Qubit[], key : Qubit[], flag : Qubit, phase : Qubit[]
    ) : Unit is Adj + Ctl {
        PrepareCandidateExperiment(db, target, pairMode, candidate, sample, key, flag);
        for q in phase { H(q); }
        for bit in 0..Length(phase)-1 {
            for _ in 1..(1 <<< bit) {
                Controlled CandidateGrover([phase[bit]], (db, target, pairMode, candidate, sample, key, flag));
            }
        }
        Adjoint FullQft(phase);
    }

    function PhaseEstimateProbability(value : Int, bits : Int) : Double {
        let angle = PI() * IntAsDouble(value) / IntAsDouble(1 <<< bits);
        let s = Sin(angle);
        return s * s;
    }

    function ContainmentToJaccard(p : Double, a : Int, b : Int) : Double {
        let upper = MinD(1.0, IntAsDouble(b) / IntAsDouble(a));
        let clipped = MaxD(0.0, MinD(upper, p));
        let intersection = IntAsDouble(a) * clipped;
        return intersection / (IntAsDouble(a + b) - intersection);
    }

    operation CandidateEstimateOnce(db : Int[][], target : Int[], pairMode : Bool,
        encoded : Int, phaseBits : Int) : (Double, Int) {
        let rbits = BitsNeeded(Length(db));
        let cbits = pairMode ? 2 * rbits | rbits;
        use candidate = Qubit[cbits];
        use sample = Qubit[DataSampleWidth(db)];
        use key = Qubit[DataKeyWidth(db, target)];
        use flag = Qubit();
        use phase = Qubit[phaseBits];
        XorConstant(encoded, candidate);
        CoherentCandidateQae(db, target, pairMode, candidate, sample, key, flag, phase);
        let outcome = MeasureInteger(phase);
        ResetAll(candidate + sample + key + [flag]);
        mutable p = PhaseEstimateProbability(outcome, phaseBits);
        if pairMode {
            let mask = (1 <<< rbits) - 1;
            let i = encoded &&& mask;
            let j = encoded >>> rbits;
            set p = ContainmentToJaccard(p, Length(db[i]), Length(db[j]));
        }
        return (p, 2 * (1 <<< phaseBits) - 1);
    }

    function Median(values : Double[]) : Double {
        mutable sorted = values;
        for i in 1..Length(sorted)-1 {
            mutable j = i;
            while j > 0 and sorted[j] < sorted[j-1] {
                let tmp = sorted[j];
                set sorted w/= j <- sorted[j-1];
                set sorted w/= j-1 <- tmp;
                set j -= 1;
            }
        }
        return sorted[Length(sorted)/2];
    }

    operation CandidateEstimate(db : Int[][], target : Int[], pairMode : Bool,
        encoded : Int, phaseBits : Int, repetitions : Int) : (Double, Int) {
        if phaseBits < 1 or phaseBits > 12 { fail "phaseBits must lie in 1..12 for these small experiments."; }
        if repetitions < 1 or repetitions % 2 == 0 { fail "Use a positive odd median repetition count."; }
        mutable estimates = [];
        mutable calls = 0;
        for _ in 1..repetitions {
            let (p, count) = CandidateEstimateOnce(db, target, pairMode, encoded, phaseBits);
            set estimates += [p];
            set calls += count;
        }
        return (Median(estimates), calls);
    }

    operation EstimateContainment(source : Int[], target : Int[], phaseBits : Int,
        repetitions : Int) : (Double, Int) {
        ValidateSet(source, false);
        ValidateSet(target, true);
        return CandidateEstimate([source], target, false, 0, phaseBits, repetitions);
    }

    operation EstimateJaccard(source : Int[], target : Int[], phaseBits : Int,
        repetitions : Int) : (Double, Int) {
        let (p, calls) = EstimateContainment(source, target, phaseBits, repetitions);
        return (ContainmentToJaccard(p, Length(source), Length(target)), calls);
    }

    // Exhaustive basis-state diagnostic of key lookup, exact membership and uncomputation.
    // This is a TEST, not a similarity-estimation algorithm.
    operation CheckMembershipOracle(source : Int[], target : Int[]) : Bool {
        ValidateSet(source, false); ValidateSet(target, true);
        use address = Qubit[BitsNeeded(Length(source))];
        use key = Qubit[DataKeyWidth([source], target)];
        use flag = Qubit();
        mutable correct = true;
        for i in 0..Length(source)-1 {
            XorConstant(i, address);
            XorLookup(source, address, key);
            Membership(target, key, flag);
            let measured = MResetZ(flag);
            mutable expected = false;
            for x in target { if x == source[i] { set expected = true; } }
            if (measured == One) != expected { set correct = false; }
            Adjoint XorLookup(source, address, key);
            XorConstant(i, address);
        }
        // Returning with all registers zero is also checked by the simulator.
        return correct;
    }
}
