using Statistics

const NSAMP = 2000
const BASELINE = 20.0
const VARIANCE = 9.0
const PROBE_OVERHEAD = 3.0
const SLO_P95 = 26.0
const SEED = 101

function run_lr11(seed::Integer)
    engee.set_param!("LR11_passive_qos/JitterNoise", "Seed" => string(seed))
    engee.run("LR11_passive_qos")
    passive = collect(Float64, passive_delay.value)
    active = collect(Float64, active_delay.value)
    p95_passive = quantile(passive, 0.95)
    p95_active = quantile(active, 0.95)
    return (passive=passive, active=active, p95_passive=p95_passive, p95_active=p95_active)
end

r1 = run_lr11(SEED)
r2 = run_lr11(SEED)
r3 = run_lr11(SEED + 1000)
r4 = run_lr11(SEED)

EQ_1_2p = r1.passive == r2.passive
EQ_1_4p = r1.passive == r4.passive
EQ_1_3p = r1.passive == r3.passive
EQ_2_4p = r2.passive == r4.passive

EQ_1_2a = r1.active == r2.active
EQ_1_4a = r1.active == r4.active
EQ_1_3a = r1.active == r3.active
EQ_2_4a = r2.active == r4.active

sigma = sqrt(VARIANCE)
p95_passive_theory = BASELINE + 1.6448536269514722 * sigma
p95_active_theory = BASELINE + PROBE_OVERHEAD + 1.6448536269514722 * sigma

bias_meas = mean(r1.active) - mean(r1.passive)

println("p95_passive=", r1.p95_passive, " p95_active=", r1.p95_active)
println("p95_passive_theory=", p95_passive_theory, " p95_active_theory=", p95_active_theory)
println("bias_meas=", bias_meas, " probe_overhead=", PROBE_OVERHEAD)
println("SLO_p95=", SLO_P95, " SLO_MET_passive=", r1.p95_passive < SLO_P95, " SLO_MET_active=", r1.p95_active < SLO_P95)
println("EQ_passive: ", EQ_1_2p, " ", EQ_1_4p, " ", EQ_1_3p, " ", EQ_2_4p)
println("EQ_active: ", EQ_1_2a, " ", EQ_1_4a, " ", EQ_1_3a, " ", EQ_2_4a)

ratio_p = r1.p95_passive / p95_passive_theory
ratio_a = r1.p95_active / p95_active_theory
PASS = EQ_1_2p && EQ_1_4p && !EQ_1_3p && EQ_2_4p &&
       EQ_1_2a && EQ_1_4a && !EQ_1_3a && EQ_2_4a &&
       0.8 <= ratio_p <= 1.2 && 0.8 <= ratio_a <= 1.2 &&
       isapprox(bias_meas, PROBE_OVERHEAD; atol=1e-9)
println("PASS=", PASS)
