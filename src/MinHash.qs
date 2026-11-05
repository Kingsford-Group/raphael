namespace QuantumSketchExperiments {
    open Std.Arithmetic;
    open Std.Canon;
    open Std.Convert;
    open Std.Measurement;
    open Std.Random;

    // This module uses an explicit multiplexed QROM circuit for functional
    // simulation. Its gate count is not the ideal-QRAM query cost in the paper.
    // Restriction: n and modulus are powers of two, n >= 2, and all keys,
    // multipliers, and offsets are nonnegative and smaller than modulus.

    function MinHashLog2PowerOfTwo(value : Int) : Int {
        mutable bits = 0;
        mutable capacity = 1;
        while capacity < value {
            set capacity *= 2;
            set bits += 1;
        }
        if capacity != value { fail "Expected a positive power of two."; }
        return bits;
    }

    function MinHashCeilingSqrt(value : Int) : Int {
        mutable root = 1;
        while root * root < value { set root += 1; }
        return root;
    }

    // A reversible, explicit lookup of input keys, not hash values.
    // |i>|z> -> |i>|z XOR keys[i]>.
    operation MinHashLoadKey(keys : Int[], index : Qubit[], key : Qubit[]) : Unit is Adj + Ctl {
        for i in 0..Length(keys) - 1 {
            for bit in 0..Length(key) - 1 {
                if ((keys[i] >>> bit) &&& 1) == 1 {
                    ApplyControlledOnInt(i, X, index, key[bit]);
                }
            }
        }
    }

    // Compute an affine hash into a zero-initialized output register using
    // reversible addition controlled by the input key bits. The arithmetic is
    // modulo 2^Length(output), so no list of hash values is precomputed.
    operation MinHashComputeHash(multiplier : Int, offset : Int, key : Qubit[], output : Qubit[]) : Unit is Adj + Ctl {
        let modulus = 1 <<< Length(output);
        IncByI(offset, output);
        for bit in 0..Length(key) - 1 {
            let increment = (multiplier * (1 <<< bit)) % modulus;
            if increment != 0 {
                Controlled IncByI([key[bit]], (increment, output));
            }
        }
    }

    operation MinHashMarkLess(keys : Int[], multiplier : Int, offset : Int,
                              modulus : Int, threshold : Int, index : Qubit[]) : Unit is Adj + Ctl {
        let width = MinHashLog2PowerOfTwo(modulus);
        use key = Qubit[width];
        use hashed = Qubit[width];
        use phase = Qubit();
        X(phase);
        within {
            MinHashLoadKey(keys, index, key);
            MinHashComputeHash(multiplier, offset, key, hashed);
        } apply {
            // A conditional Z on |1> implements the required phase flip.
            ApplyIfGreaterL(Z, IntAsBigInt(threshold), hashed, phase);
        }
        X(phase);
    }

    operation MinHashDiffuse(index : Qubit[]) : Unit is Adj + Ctl {
        for q in index { H(q); X(q); }
        let last = Length(index) - 1;
        if last == 0 {
            Z(index[0]);
        } else {
            Controlled Z(index[0..last - 1], index[last]);
        }
        for q in index { X(q); H(q); }
    }

    operation MinHashSampleCandidate(keys : Int[], multiplier : Int, offset : Int,
                                     modulus : Int, threshold : Int, iterations : Int) : Int {
        use index = Qubit[MinHashLog2PowerOfTwo(Length(keys))];
        for q in index { H(q); }
        for _ in 1..iterations {
            MinHashMarkLess(keys, multiplier, offset, modulus, threshold, index);
            MinHashDiffuse(index);
        }
        mutable candidate = 0;
        for bit in 0..Length(index) - 1 {
            if MResetZ(index[bit]) == One { set candidate += 1 <<< bit; }
        }
        return candidate;
    }

    // Evaluate the measured candidate through the same quantum data/hash
    // circuit. No classical exhaustive minimum is used by the search.
    operation MinHashReadValue(keys : Int[], multiplier : Int, offset : Int,
                               modulus : Int, candidate : Int) : Int {
        let width = MinHashLog2PowerOfTwo(modulus);
        use index = Qubit[MinHashLog2PowerOfTwo(Length(keys))];
        use key = Qubit[width];
        use hashed = Qubit[width];
        within {
            ApplyXorInPlace(candidate, index);
            MinHashLoadKey(keys, index, key);
        } apply {
            MinHashComputeHash(multiplier, offset, key, hashed);
        }
        mutable value = 0;
        for bit in 0..width - 1 {
            if MResetZ(hashed[bit]) == One { set value += 1 <<< bit; }
        }
        return value;
    }

    // One bounded Durr-Hoyer run with a BBHT-style expanding random-iteration
    // schedule. The stopping budget is 24*ceil(sqrt(n)) Grover/evaluation units.
    // The second result counts logical hash evaluations, including uncomputation:
    // two per phase-oracle application and one per measured value evaluation.
    operation MinHashSingleRun(keys : Int[], multiplier : Int, offset : Int,
                              modulus : Int) : (Int, Int) {
        let n = Length(keys);
        let limit = MinHashCeilingSqrt(n);
        let budget = 24 * limit;
        let initial = DrawRandomInt(0, n - 1);
        mutable best = MinHashReadValue(keys, multiplier, offset, modulus, initial);
        mutable evaluations = 1;
        mutable spent = 0;
        mutable window = 1;
        while spent < budget and best > 0 {
            mutable iterations = DrawRandomInt(0, window - 1);
            if iterations + 1 > budget - spent {
                set iterations = budget - spent - 1;
            }
            let candidate = MinHashSampleCandidate(keys, multiplier, offset, modulus, best, iterations);
            let value = MinHashReadValue(keys, multiplier, offset, modulus, candidate);
            set evaluations += 2 * iterations + 1;
            set spent += iterations + 1;
            if value < best {
                set best = value;
                set window = 1;
            } else {
                // ceil(6*window/5), capped at ceil(sqrt(n)).
                set window = (6 * window + 4) / 5;
                if window > limit { set window = limit; }
            }
        }
        return (best, evaluations);
    }

    /// Construct all prescribed affine-hash minima by independent quantum
    /// minimum-finding restarts, reusing the same classical input description.
    /// Returns (signature, logical hash evaluations including inverses).
    /// The simulator/runner must state the empirical failure rate and seeds;
    /// repetitions here is a user-specified restart count, not a claimed delta.
    operation MultiHashMinima(keys : Int[], multipliers : Int[], offsets : Int[],
                              modulus : Int, repetitions : Int) : (Int[], Int) {
        if Length(keys) < 2 { fail "At least two keys are required."; }
        let indexBits = MinHashLog2PowerOfTwo(Length(keys));
        let hashBits = MinHashLog2PowerOfTwo(modulus);
        if modulus < 2 or repetitions < 1 { fail "modulus >= 2 and repetitions >= 1 required."; }
        if Length(multipliers) != Length(offsets) or Length(multipliers) == 0 {
            fail "Nonempty multiplier and offset arrays must have equal lengths.";
        }
        for key in keys {
            if key < 0 or key >= modulus { fail "A key is outside the hash word range."; }
        }
        mutable signature = [];
        mutable evaluations = 0;
        for r in 0..Length(multipliers) - 1 {
            let multiplier = multipliers[r];
            let offset = offsets[r];
            if multiplier < 0 or multiplier >= modulus or offset < 0 or offset >= modulus {
                fail "A hash seed is outside the hash word range.";
            }
            mutable best = modulus;
            for _ in 1..repetitions {
                let (candidate, calls) = MinHashSingleRun(keys, multiplier, offset, modulus);
                if candidate < best { set best = candidate; }
                set evaluations += calls;
            }
            set signature += [best];
        }
        return (signature, evaluations);
    }
}
