using Statistics
using SpecialFunctions: erfc

const NSAMP = 2000
const Z_LOSS = 1.645
const SEED = 101

Qfun(x) = 0.5 * erfc(x / sqrt(2))

function run_lr12(seed::Integer)
    engee.set_param!("LR12_streaming_telemetry/MeasNoise", "Seed" => string(seed))
    engee.set_param!("LR12_streaming_telemetry/LossNoise", "Seed" => string(seed))
    engee.run("LR12_streaming_telemetry")
    raw = collect(Float64, telemetry_raw.value)
    lossg = collect(Float64, loss_gaussian.value)
    loss_flag = lossg .> Z_LOSS
    n_lost = count(loss_flag)
    completeness = 1 - n_lost / NSAMP
    max_gap = 0
    cur_gap = 0
    for f in loss_flag
        if f
            cur_gap += 1
            max_gap = max(max_gap, cur_gap)
        else
            cur_gap = 0
        end
    end
    return (raw=raw, lossg=lossg, n_lost=n_lost, completeness=completeness, max_gap=max_gap)
end

r1 = run_lr12(SEED)
r2 = run_lr12(SEED)
r3 = run_lr12(SEED + 1000)
r4 = run_lr12(SEED)

EQ_raw_1_2 = r1.raw == r2.raw
EQ_raw_1_4 = r1.raw == r4.raw
EQ_raw_1_3 = r1.raw == r3.raw
EQ_raw_2_4 = r2.raw == r4.raw

EQ_loss_1_2 = r1.lossg == r2.lossg
EQ_loss_1_4 = r1.lossg == r4.lossg
EQ_loss_1_3 = r1.lossg == r3.lossg
EQ_loss_2_4 = r2.lossg == r4.lossg

p_loss_theory = Qfun(Z_LOSS)
completeness_theory = 1 - p_loss_theory

println("n_lost=", r1.n_lost, " completeness=", r1.completeness, " completeness_theory=", completeness_theory)
println("max_gap=", r1.max_gap)
println("EQ_raw: ", EQ_raw_1_2, " ", EQ_raw_1_4, " ", EQ_raw_1_3, " ", EQ_raw_2_4)
println("EQ_loss: ", EQ_loss_1_2, " ", EQ_loss_1_4, " ", EQ_loss_1_3, " ", EQ_loss_2_4)

ratio = (1 - r1.completeness) / p_loss_theory
PASS = EQ_raw_1_2 && EQ_raw_1_4 && !EQ_raw_1_3 && EQ_raw_2_4 &&
       EQ_loss_1_2 && EQ_loss_1_4 && !EQ_loss_1_3 && EQ_loss_2_4 &&
       0.5 <= ratio <= 2.0
println("PASS=", PASS)
