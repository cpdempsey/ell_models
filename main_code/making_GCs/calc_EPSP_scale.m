function GC_model = calc_EPSP_scale(GC_model)


% tau_s = GC_model.tau_s;
% tau_m = GC_model.tau_m;
% tau_s_slow = GC_model.tau_s_slow;
tran=GC_model.dt:GC_model.dt:200;

kernel_fast = 1/(GC_model.tau_s - GC_model.tau_m)*(exp(-tran/GC_model.tau_s) - exp(-tran/GC_model.tau_m));
kernel_slow = 1/(GC_model.tau_s_slow - GC_model.tau_m)*(exp(-tran/GC_model.tau_s_slow) - exp(-tran/GC_model.tau_m));

GC_model.Wscale_fast = 1/max(kernel_fast);
GC_model.Wscale_slow = 1/max(kernel_slow);

