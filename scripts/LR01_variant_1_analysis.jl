# ЛР-01. Обработка результатов измерительного эксперимента
# Вариант 1
#
# Перед запуском скрипта выполните моделирование файла
# LR01_student_template_variant_1_final.engee. В рабочей области должны
# появиться переменные measured_data и true_data.

using Statistics
using Plots

# Параметры варианта 1
X0 = 10.0            # истинное значение, мВ
delta_sys = 0.20     # систематическая погрешность, мВ
sigma = 0.50         # заданное СКО случайной погрешности, мВ
sample_time = 0.1    # период дискретизации, с

# Преобразование WorkspaceArray в обычные массивы
if !@isdefined(measured_data)
    error("Переменная measured_data не найдена. Сначала запустите модель Engee.")
end

measured_table = collect(measured_data)
x = Float64.(collect(measured_table.value))

if length(x) < 2
    error("Для расчёта выборочного СКО требуется не менее двух отсчётов.")
end

# Временная шкала. Если WorkspaceArray содержит поле time, используем его;
# иначе формируем шкалу по заданному периоду дискретизации.
if :time in propertynames(measured_table)
    t = Float64.(collect(measured_table.time))
else
    t = collect(0:length(x)-1) .* sample_time
end

# Истинный сигнал для совместной временной диаграммы
if @isdefined(true_data)
    true_table = collect(true_data)
    x_true = Float64.(collect(true_table.value))
else
    x_true = fill(X0, length(x))
end

if length(x_true) != length(x)
    x_true = fill(X0, length(x))
end

# Статистические показатели
n = length(x)
x_mean = mean(x)
s = std(x)                    # выборочное СКО: деление на n - 1
u_A = s / sqrt(n)

# Для базовой серии варианта 1: n = 100, число степеней свободы = 99
t_crit = 1.984
ci_low = x_mean - t_crit * u_A
ci_high = x_mean + t_crit * u_A
ci_width = ci_high - ci_low

bias = x_mean - X0
relative_error = bias / X0 * 100
x_corrected = x_mean - delta_sys
residual_error = x_corrected - X0

# Контрольные критерии из методики
s_relative_deviation = abs(s - sigma) / sigma
mean_deviation = abs(x_mean - (X0 + delta_sys))
mean_limit = 3 * sigma / sqrt(n)
check_s = s_relative_deviation <= 0.20
check_mean = mean_deviation <= mean_limit

println("ЛР-01. Вариант 1")
println("Количество отсчётов n = ", n)
println("Среднее значение = ", x_mean, " мВ")
println("Выборочное СКО s = ", s, " мВ")
println("Стандартная неопределённость u_A = ", u_A, " мВ")
println("95%-й доверительный интервал = [", ci_low, "; ", ci_high, "] мВ")
println("Ширина доверительного интервала = ", ci_width, " мВ")
println("Смещение относительно X0 = ", bias, " мВ")
println("Относительная погрешность = ", relative_error, " %")
println("Результат после поправки = ", x_corrected, " мВ")
println("Остаточная погрешность после поправки = ", residual_error, " мВ")
println("Проверка СКО (отклонение не более 20 %) = ", check_s)
println("Проверка среднего по границе 3σ/√n = ", check_mean)

# Рисунок 1. Временная диаграмма истинного и измеренного сигналов
p_time = plot(
    t,
    x,
    label = "x(t) — измеренный сигнал",
    xlabel = "t, с",
    ylabel = "Напряжение, мВ",
    title = "ЛР-01, вариант 1: временная диаграмма",
    linewidth = 2,
    grid = true,
    legend = :topright,
)
plot!(p_time, t, x_true, label = "X0 — истинное значение", linewidth = 2)
display(p_time)

# Рисунок 2. Гистограмма базовой серии
p_hist = histogram(
    x,
    bins = 10,
    label = "Результаты измерений",
    xlabel = "x, мВ",
    ylabel = "Частота",
    title = "ЛР-01, вариант 1: гистограмма результатов",
    grid = true,
    legend = :topright,
)
vline!(p_hist, [x_mean], label = "Среднее", linewidth = 2)
vline!(p_hist, [X0], label = "X0", linewidth = 2)
display(p_hist)

