using Statistics

const NSAMP = 2000
const TRAIN_N = 1000
const ANOMALY_START = 10.0
const DT = 0.01
const K_SIGMA = 3.0
const SEED = 101

function run_lr16(seed)
    engee.set_param!("LR16_cross_channel_diagnostics/DelayJitter", "Seed" => string(seed))
    engee.set_param!("LR16_cross_channel_diagnostics/MeasNoise", "Seed" => string(seed))
    engee.set_param!("LR16_cross_channel_diagnostics/ControlNoise", "Seed" => string(seed + 300))
    engee.run("LR16_cross_channel_diagnostics")
    delay = collect(Float64, delay_metric.value)
    telemetry = collect(Float64, telemetry_value.value)
    control = collect(Float64, control_metric.value)
    return (delay=delay, telemetry=telemetry, control=control)
end

function threshold_flag(x)
    train = x[1:TRAIN_N]
    bmean = mean(train)
    bstd = std(train)
    thr = bmean + K_SIGMA * bstd
    return (flag=(x .> thr), bmean=bmean, bstd=bstd, thr=thr)
end

r1 = run_lr16(SEED)
r2 = run_lr16(SEED)
r3 = run_lr16(SEED + 1000)
r4 = run_lr16(SEED)

is_anomaly_period = [((k - 1) * DT) >= ANOMALY_START for k in 1:NSAMP]
n_anomaly = count(is_anomaly_period)
n_normal = NSAMP - n_anomaly

fd = threshold_flag(r1.delay)
ft = threshold_flag(r1.telemetry)
fc = threshold_flag(r1.control)

tpr_d = count(fd.flag .& is_anomaly_period) / n_anomaly
fpr_d = count(fd.flag .& .!is_anomaly_period) / n_normal
tpr_t = count(ft.flag .& is_anomaly_period) / n_anomaly
fpr_t = count(ft.flag .& .!is_anomaly_period) / n_normal
tpr_c = count(fc.flag .& is_anomaly_period) / n_anomaly
fpr_c = count(fc.flag .& .!is_anomaly_period) / n_normal

joint_flag = fd.flag .& ft.flag
tpr_joint = count(joint_flag .& is_anomaly_period) / n_anomaly
fpr_joint = count(joint_flag .& .!is_anomaly_period) / n_normal

corr_dt = cor(r1.delay, r1.telemetry)
corr_dc = cor(r1.delay, r1.control)

EQ_delay_1_2 = r1.delay == r2.delay
EQ_delay_1_4 = r1.delay == r4.delay
EQ_delay_1_3 = r1.delay == r3.delay
EQ_delay_2_4 = r2.delay == r4.delay

EQ_tel_1_2 = r1.telemetry == r2.telemetry
EQ_tel_1_3 = r1.telemetry == r3.telemetry

EQ_ctrl_1_2 = r1.control == r2.control
EQ_ctrl_1_3 = r1.control == r3.control

println("N_delay=", length(r1.delay), " N_telemetry=", length(r1.telemetry), " N_control=", length(r1.control))
println("tpr_delay=", tpr_d, " fpr_delay=", fpr_d)
println("tpr_telemetry=", tpr_t, " fpr_telemetry=", fpr_t)
println("tpr_control=", tpr_c, " fpr_control=", fpr_c)
println("tpr_joint=", tpr_joint, " fpr_joint=", fpr_joint)
println("corr_delay_telemetry=", corr_dt, " corr_delay_control=", corr_dc)
println("EQ_delay_1_2=", EQ_delay_1_2, " EQ_delay_1_4=", EQ_delay_1_4, " EQ_delay_1_3=", EQ_delay_1_3, " EQ_delay_2_4=", EQ_delay_2_4)
println("EQ_tel_1_2=", EQ_tel_1_2, " EQ_tel_1_3=", EQ_tel_1_3)
println("EQ_ctrl_1_2=", EQ_ctrl_1_2, " EQ_ctrl_1_3=", EQ_ctrl_1_3)

PASS = length(r1.delay) == NSAMP && length(r1.telemetry) == NSAMP && length(r1.control) == NSAMP &&
       EQ_delay_1_2 && EQ_delay_1_4 && !EQ_delay_1_3 && EQ_delay_2_4 &&
       EQ_tel_1_2 && !EQ_tel_1_3 &&
       EQ_ctrl_1_2 && !EQ_ctrl_1_3 &&
       tpr_joint > 0.7 && fpr_joint < 0.05 &&
       tpr_c < 0.1 &&
       corr_dt > 0.3 && abs(corr_dc) < 0.3
println("PASS=", PASS)
