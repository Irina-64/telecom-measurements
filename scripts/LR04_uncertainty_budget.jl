using Statistics

const MODEL_NAME = "LR04_uncertainty_budget"
const RANDOM_BLOCK_PATH = "LR04_uncertainty_budget/RandomNoise"
const BASE_SEED = 101
const ALT_SEED = BASE_SEED + 1000
const EXPECTED_LENGTH = 100
const X0 = 10.0
const U_REF = 0.05
const K_COVERAGE = 2.0

function run_model_and_extract()
    engee.run(MODEL_NAME)
        x = Float64.(collect(measured_data).value)
            @assert !isempty(x) "Выход measured_data пуст"
    @assert length(x) == EXPECTED_LENGTH "Неожиданная длина выхода"
        @assert all(isfinite, x) "Выход содержит NaN или Inf"
    return x
    end

function set_seed!(seed::Integer)
    engee.set_param!(RANDOM_BLOCK_PATH, "Seed" => string(seed))
end

function budget_metrics(x)
    n = length(x)
        xbar = mean(x)
            s = std(x)
                u_a = s / sqrt(n)
    u_b = U_REF / sqrt(3)
    u_c = sqrt(u_a^2 + u_b^2)
        U_exp = K_COVERAGE * u_c
    bias_hat = xbar - X0
        correction = -bias_hat
            x_corr = xbar + correction
                covers_x0 = abs(x_corr - X0) <= U_exp
                    return (n=n, xbar=xbar, s=s, u_a=u_a, u_b=u_b, u_c=u_c, U_exp=U_exp, bias_hat=bias_hat, correction=correction, x_corr=x_corr, covers_x0=covers_x0)
                    end

function check_numeric_acceptance(m)
    @assert isapprox(m.bias_hat, 0.20; atol=0.15) "Смещение вне допуска"
    @assert 0.3 <= m.s <= 0.7 "СКО вне допуска"
        @assert m.covers_x0 "Опорное значение не попадает в охватывающий интервал"
            return true
            end

println("MODEL = ", MODEL_NAME)
println("SEEDS = ", (BASE_SEED, BASE_SEED, ALT_SEED, BASE_SEED))

set_seed!(BASE_SEED)
x1 = run_model_and_extract()
m1 = budget_metrics(x1)

set_seed!(BASE_SEED)
x2 = run_model_and_extract()
m2 = budget_metrics(x2)

set_seed!(ALT_SEED)
x3 = run_model_and_extract()
m3 = budget_metrics(x3)

set_seed!(BASE_SEED)
x4 = run_model_and_extract()
m4 = budget_metrics(x4)

eq_1_2 = x1 == x2
eq_1_4 = x1 == x4
eq_1_3 = x1 == x3
eq_2_4 = x2 == x4

@assert eq_1_2 "Запуски 1 и 2 при одинаковом seed различаются"
@assert eq_1_4 "После восстановления seed выход не воспроизведён"
@assert !eq_1_3 "Изменение seed не изменило полный выходной массив"
@assert eq_2_4 "Запуски 2 и 4 при одинаковом seed различаются"

check_numeric_acceptance(m1)
check_numeric_acceptance(m2)
check_numeric_acceptance(m4)

set_seed!(BASE_SEED)
engee.save(MODEL_NAME, MODEL_NAME * ".engee"; force=true)

println("RUN_1 = ", m1)
println("RUN_2 = ", m2)
println("RUN_3 = ", m3)
println("RUN_4 = ", m4)
println("EQ_1_2 = ", eq_1_2)
println("EQ_1_4 = ", eq_1_4)
println("EQ_1_3 = ", eq_1_3)
println("EQ_2_4 = ", eq_2_4)
println("PASS: модель и критерии варианта 1 проверены; параметры восстановлены")
