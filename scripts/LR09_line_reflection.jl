using Statistics

const NSAMP = 2000
const DT = 0.001
const A0 = 2.0
const PERIOD = 0.2
const PULSE_WIDTH_PCT = 15
const DELAY = 0.05
const GAMMA_THEORY = 0.5
const SEED = 101

function reflected_window_mean(x, dt::Float64, period::Float64, duty_pct::Real, delay::Float64, n::Int)
    dur = (duty_pct / 100) * period
    total = 0.0
    count = 0
    for m in 0:n-1
        t = m * dt
        tau = mod(t, period)
        if delay <= tau < delay + dur
            total += x[m+1]
            count += 1
        end
    end
    return total / count
end

function run_lr09(seed::Integer)
    engee.set_param!("LR09_line_reflection/MeasNoise", "Seed" => string(seed))
    engee.run("LR09_line_reflection")
    data = copy(measured_reflection.value)
    a_reflected = reflected_window_mean(data, DT, PERIOD, PULSE_WIDTH_PCT, DELAY, NSAMP)
    return (data=data, a_reflected=a_reflected)
end

r1 = run_lr09(SEED)
r2 = run_lr09(SEED)
r3 = run_lr09(SEED + 1000)
r4 = run_lr09(SEED)

EQ_1_2 = r1.data == r2.data
EQ_1_4 = r1.data == r4.data
EQ_1_3 = r1.data == r3.data
EQ_2_4 = r2.data == r4.data

gamma_meas = r1.a_reflected / A0
RL_meas = -20 * log10(abs(gamma_meas))
RL_theory = -20 * log10(abs(GAMMA_THEORY))

println("a_reflected=", r1.a_reflected, " gamma_meas=", gamma_meas, " gamma_theory=", GAMMA_THEORY)
println("RL_meas=", RL_meas, " RL_theory=", RL_theory)
println("EQ_1_2=", EQ_1_2, " EQ_1_4=", EQ_1_4, " EQ_1_3=", EQ_1_3, " EQ_2_4=", EQ_2_4)

ratio = gamma_meas / GAMMA_THEORY
PASS = EQ_1_2 && EQ_1_4 && !EQ_1_3 && EQ_2_4 && 0.5 <= ratio <= 2.0
println("PASS=", PASS)
