using Statistics
using SpecialFunctions: erfc

const A0 = 2.0
const PERIOD = 0.2
const PULSE_WIDTH_PCT = 30
const PHASE_DELAY = 0.0
const DT = 0.001
const NSAMP = 2000
const THRESHOLD = A0 / 2
const VARIANCE = 0.16
const SEED = 101

Qfun(x) = 0.5 * erfc(x / sqrt(2))

function ideal_high(t, T, duty_pct, phase)
    tau = mod(t - phase, T)
    return tau < (duty_pct / 100) * T
end

function run_lr06(seed::Integer)
    engee.set_param!("LR06_threshold_ber/RandomNoise", "Seed" => string(seed))
    engee.run("LR06_threshold_ber")
    data = copy(measured_data.value)
    t = collect(0:DT:(NSAMP-1)*DT)
    ideal = [ideal_high(ti, PERIOD, PULSE_WIDTH_PCT, PHASE_DELAY) ? A0 : 0.0 for ti in t]
    decision = data .> THRESHOLD
    truth = ideal .> THRESHOLD
    errors = count(decision .!= truth)
    ber_hat = errors / NSAMP
    return (data=data, errors=errors, ber_hat=ber_hat)
end

r1 = run_lr06(SEED)
r2 = run_lr06(SEED)
r3 = run_lr06(SEED + 1000)
r4 = run_lr06(SEED)

EQ_1_2 = r1.data == r2.data
EQ_1_4 = r1.data == r4.data
EQ_1_3 = r1.data == r3.data
EQ_2_4 = r2.data == r4.data

d = THRESHOLD
sigma = sqrt(VARIANCE)
ber_theory = Qfun(d / sigma)

ber_metrics = (
    N = NSAMP,
    errors = r1.errors,
    ber_hat = r1.ber_hat,
    ber_theory = ber_theory,
    d = d,
    sigma = sigma,
    d_over_sigma = d / sigma,
)

println("errors=", r1.errors, " ber_hat=", r1.ber_hat, " ber_theory=", ber_theory)
println("EQ_1_2=", EQ_1_2, " EQ_1_4=", EQ_1_4, " EQ_1_3=", EQ_1_3, " EQ_2_4=", EQ_2_4)

PASS = EQ_1_2 && EQ_1_4 && !EQ_1_3 && EQ_2_4 && r1.errors > 0
println("PASS=", PASS)
