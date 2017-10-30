function vol = volume(metrics,op)
  vol = metrics.J*kron(op.Px,op.Py);
end
