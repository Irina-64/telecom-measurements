using Statistics

const NSAMP = 2000
const SEED = 101
const Z_THRESH = 1.645

function run_lr15(seed)
    engee.set_param!("LR15_integrative_measurement/DelayJitter", "Seed" => string(seed))
    engee.set_param!("LR15_integrative_measurement/MeasNoise", "Seed" => string(seed))
    engee.set_param!("LR15_integrative_measurement/LossNoise", "Seed" => string(seed))
    engee.run("LR15_integrative_measurement")
    delay = collect(Float64, delay_metric.value)
    telemetry = collect(Float64, telemetry_value.value)
    loss = collect(Float64, loss_gaussian.value)
    return (delay=delay, telemetry=telemetry, loss=loss)
end

r1 = run_lr15(SEED)
r2 = run_lr15(SEED)
r3 = run_lr15(SEED + 1000)
r4 = run_lr15(SEED)

EQ_delay_1_2 = r1.delay == r2.delay
EQ_delay_1_4 = r1.delay == r4.delay
EQ_delay_1_3 = r1.delay == r3.delay
EQ_delay_2_4 = r2.delay == r4.delay

EQ_tel_1_2 = r1.telemetry == r2.telemetry
EQ_tel_1_3 = r1.telemetry == r3.telemetry

EQ_loss_1_2 = r1.loss == r2.loss
EQ_loss_1_3 = r1.loss == r3.loss

p50 = quantile(r1.delay, 0.50)
p95 = quantile(r1.delay, 0.95)

n_lost = count(x -> x > Z_THRESH, r1.loss)
completeness = 1 - n_lost / NSAMP

println("N_delay=", length(r1.delay), " N_telemetry=", length(r1.telemetry), " N_loss=", length(r1.loss))
println("p50=", p50, " p95=", p95)
println("n_lost=", n_lost, " completeness=", completeness)
println("EQ_delay_1_2=", EQ_delay_1_2, " EQ_delay_1_4=", EQ_delay_1_4, " EQ_delay_1_3=", EQ_delay_1_3, " EQ_delay_2_4=", EQ_delay_2_4)
println("EQ_tel_1_2=", EQ_tel_1_2, " EQ_tel_1_3=", EQ_tel_1_3)
println("EQ_loss_1_2=", EQ_loss_1_2, " EQ_loss_1_3=", EQ_loss_1_3)

PASS = length(r1.delay) == NSAMP && length(r1.telemetry) == NSAMP && length(r1.loss) == NSAMP &&
       EQ_delay_1_2 && EQ_delay_1_4 && !EQ_delay_1_3 && EQ_delay_2_4 &&
       EQ_tel_1_2 && !EQ_tel_1_3 &&
       EQ_loss_1_2 && !EQ_loss_1_3 &&
       abs(p95 - 24.935) < 3.0 && completeness > 0.9
println("PASS=", PASS)
