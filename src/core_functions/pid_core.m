function pid_out = pid_core(pdf, opts)
pid_out = feval(opts.function, pdf);
end