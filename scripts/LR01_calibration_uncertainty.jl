using Statistics

const NSAMP = 100
const SEED = 101
const X0 = 10.0
const DELTA_SYS = 0.20
const SIGMA = 0.5

function run_lr01(seed)
    engee.set_param!("LR01_calibration_uncertainty/RandomNoise", "Seed" => string(seed))
    engee.run("LR01_calibration_uncertainty")
    return collect(Float64, measured_data.value)
end

r1 = run_lr01(SEED)
r2 = run_lr01(SEED)
r3 = run_lr01(SEED + 1000)
r4 = run_lr01(SEED)

EQ_1_2 = r1 == r2
EQ_1_4 = r1 == r4
EQ_1_3 = r1 == r3
EQ_2_4 = r2 == r4

n = length(r1)
x_mean = mean(r1)
s = std(r1)
u_A = s / sqrt(n)
t_crit = 1.984
ci_low = x_mean - t_crit * u_A
ci_high = x_mean + t_crit * u_A
bias = x_mean - X0
x_corrected = x_mean - DELTA_SYS
residual = x_corrected - X0

check_s = abs(s - SIGMA) / SIGMA <= 0.20
mean_limit = 3 * SIGMA / sqrt(n)
check_mean = abs(x_mean - (X0 + DELTA_SYS)) <= mean_limit

println("n=", n, " mean=", x_mean, " s=", s, " u_A=", u_A)
println("CI=[", ci_low, ", ", ci_high, "]")
println("bias=", bias, " corrected=", x_corrected, " residual=", residual)
println("EQ_1_2=", EQ_1_2, " EQ_1_4=", EQ_1_4, " EQ_1_3=", EQ_1_3, " EQ_2_4=", EQ_2_4)
println("check_s=", check_s, " check_mean=", check_mean)

PASS = n == NSAMP && EQ_1_2 && EQ_1_4 && !EQ_1_3 && EQ_2_4 && check_s && check_mean &&
       all(isfinite, (x_mean, s, u_A, bias, x_corrected, residual))
println("PASS=", PASS)
