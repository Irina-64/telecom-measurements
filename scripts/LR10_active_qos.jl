using Statistics

const NSAMP = 2000
const BASELINE = 20.0
const VARIANCE = 9.0
const SLO_P95 = 30.0
const SEED = 101

function run_lr10(seed::Integer)
    engee.set_param!("LR10_active_qos/JitterNoise", "Seed" => string(seed))
    engee.run("LR10_active_qos")
    data = collect(Float64, measured_delay.value)
    p50 = quantile(data, 0.50)
    p95 = quantile(data, 0.95)
    p99 = quantile(data, 0.99)
    return (data=data, p50=p50, p95=p95, p99=p99)
end

r1 = run_lr10(SEED)
r2 = run_lr10(SEED)
r3 = run_lr10(SEED + 1000)
r4 = run_lr10(SEED)

EQ_1_2 = r1.data == r2.data
EQ_1_4 = r1.data == r4.data
EQ_1_3 = r1.data == r3.data
EQ_2_4 = r2.data == r4.data

sigma = sqrt(VARIANCE)
p50_theory = BASELINE
p95_theory = BASELINE + 1.6448536269514722 * sigma
p99_theory = BASELINE + 2.3263478740408408 * sigma

println("p50=", r1.p50, " p95=", r1.p95, " p99=", r1.p99)
println("p50_theory=", p50_theory, " p95_theory=", p95_theory, " p99_theory=", p99_theory)
println("SLO_p95=", SLO_P95, " SLO_MET=", r1.p95 < SLO_P95)
println("EQ_1_2=", EQ_1_2, " EQ_1_4=", EQ_1_4, " EQ_1_3=", EQ_1_3, " EQ_2_4=", EQ_2_4)

ratio95 = r1.p95 / p95_theory
PASS = EQ_1_2 && EQ_1_4 && !EQ_1_3 && EQ_2_4 && 0.8 <= ratio95 <= 1.2
println("PASS=", PASS)
