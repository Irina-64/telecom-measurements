# ЛР-02. Измерение амплитуды, частоты и фазового сдвига, вариант 1
# Перед запуском выполните моделирование. Модель должна создать
# WorkspaceArray-переменные sig1_data и sig2_data.

using Statistics
using LinearAlgebra
using Printf

const A1_NOM = 1.0
const A2_NOM = 0.8
const F_NOM = 50.0
const DPHI_NOM_DEG = 30.0
const TS_NOM = 0.0002

"Преобразовать WorkspaceArray Engee в два вектора: время и значение."
function read_signal(data, signal_name::AbstractString)
    frame = collect(data)
    hasproperty(frame, :time) || error("$signal_name: отсутствует столбец time")
    hasproperty(frame, :value) || error("$signal_name: отсутствует столбец value")
    t = Float64.(frame.time)
    x = Float64.(frame.value)
    length(t) == length(x) || error("$signal_name: длины time и value различаются")
    length(t) >= 3 || error("$signal_name: недостаточно отсчётов")
    all(isfinite, t) && all(isfinite, x) || error("$signal_name: есть NaN/Inf")
    all(diff(t) .> 0) || error("$signal_name: время должно строго возрастать")
    return t, x
end

amplitude_pp(x) = (maximum(x) - minimum(x)) / 2
rms_value(x) = sqrt(mean(abs2, x))
wrap180(phi_deg) = atan(sinpi(phi_deg / 180), cospi(phi_deg / 180)) * 180 / pi

"Моменты положительных переходов через ноль с линейной интерполяцией."
function rising_zero_crossings(t, x)
    idx = findall(i -> x[i] <= 0 && x[i + 1] > 0, 1:length(x)-1)
    length(idx) >= 2 || error("Недостаточно положительных переходов через ноль")
    return [t[i] - x[i] * (t[i + 1] - t[i]) / (x[i + 1] - x[i]) for i in idx]
end

function frequency_zc(t, x)
    crossings = rising_zero_crossings(t, x)
    return 1 / mean(diff(crossings)), crossings
end

"Фаза синусоиды по МНК: x=a*sin(wt)+b*cos(wt)+c."
function phase_ls(t, x, f)
    omega = 2pi * f
    design = hcat(sin.(omega .* t), cos.(omega .* t), ones(length(t)))
    coeff = design \ x
    phase_deg = atan(coeff[2], coeff[1]) * 180 / pi
    amplitude = hypot(coeff[1], coeff[2])
    offset = coeff[3]
    return phase_deg, amplitude, offset
end

t1, x1 = read_signal(sig1_data, "sig1_data")
t2, x2 = read_signal(sig2_data, "sig2_data")
length(t1) == length(t2) || error("Число отсчётов каналов различается")
maximum(abs.(t1 .- t2)) <= 10eps(maximum(abs, t1)) || error("Шкалы времени каналов не совпадают")

dt = diff(t1)
ts_meas = mean(dt)
maximum(abs.(dt .- ts_meas)) <= 100eps(ts_meas) || @warn "Шаг времени не вполне равномерен"

f1, zc1 = frequency_zc(t1, x1)
f2, zc2 = frequency_zc(t2, x2)
phi1, als1, c1 = phase_ls(t1, x1, F_NOM)
phi2, als2, c2 = phase_ls(t2, x2, F_NOM)
dphi_ls = wrap180(phi2 - phi1)
dphi_zc = wrap180(-360F_NOM * (zc2[1] - zc1[1]))

results = (
    N = length(t1),
    t_start = first(t1),
    t_stop = last(t1),
    Ts = ts_meas,
    A1_pp = amplitude_pp(x1),
    A2_pp = amplitude_pp(x2),
    A1_ls = als1,
    A2_ls = als2,
    RMS1 = rms_value(x1),
    RMS2 = rms_value(x2),
    f1 = f1,
    f2 = f2,
    phi1_deg = wrap180(phi1),
    phi2_deg = wrap180(phi2),
    dphi_ls_deg = dphi_ls,
    dphi_zc_deg = dphi_zc,
    dc1 = c1,
    dc2 = c2,
)

@printf("ЛР-02, вариант 1\n")
@printf("N = %d, t = %.7f…%.7f с, Ts = %.7f с\n", results.N, results.t_start, results.t_stop, results.Ts)
@printf("x1: A(pp)=%.9f В, A(МНК)=%.9f В, RMS=%.9f В, f=%.9f Гц, phi=%.9f°\n",
        results.A1_pp, results.A1_ls, results.RMS1, results.f1, results.phi1_deg)
@printf("x2: A(pp)=%.9f В, A(МНК)=%.9f В, RMS=%.9f В, f=%.9f Гц, phi=%.9f°\n",
        results.A2_pp, results.A2_ls, results.RMS2, results.f2, results.phi2_deg)
@printf("Δphi: МНК=%.9f°, нулевые переходы=%.9f°, расхождение=%.9f°\n",
        results.dphi_ls_deg, results.dphi_zc_deg,
        abs(wrap180(results.dphi_ls_deg - results.dphi_zc_deg)))

# Контрольные допуски из задания.
checks = (
    amplitude_1 = abs(results.A1_pp - A1_NOM) / A1_NOM <= 0.005,
    amplitude_2 = abs(results.A2_pp - A2_NOM) / A2_NOM <= 0.005,
    frequency_1 = abs(results.f1 - F_NOM) / F_NOM <= 0.002,
    frequency_2 = abs(results.f2 - F_NOM) / F_NOM <= 0.002,
    phase_ls = abs(wrap180(results.dphi_ls_deg - DPHI_NOM_DEG)) <= 0.5,
    phase_methods = abs(wrap180(results.dphi_ls_deg - results.dphi_zc_deg)) <= 1.0,
    rms_1 = abs(results.RMS1 / A1_NOM - inv(sqrt(2))) <= 0.01 / sqrt(2),
    rms_2 = abs(results.RMS2 / A2_NOM - inv(sqrt(2))) <= 0.01 / sqrt(2),
)

println("Контроль: ", checks)
all(values(checks)) || error("Один или несколько контрольных допусков не выполнены")
println("Итог: PASS")

# Необязательная визуализация первых трёх периодов (если установлен Plots).
try
    @eval using Plots
    mask = t1 .<= first(t1) + 3 / F_NOM
    p = plot(t1[mask], x1[mask], label="x1", lw=2,
             xlabel="Время, с", ylabel="Амплитуда, В",
             title="ЛР-02, вариант 1: первые 3 периода", grid=true)
    plot!(p, t2[mask], x2[mask], label="x2", lw=2)
    display(p)
catch err
    @info "График не построен; численные расчёты завершены" exception=(err, catch_backtrace())
end

