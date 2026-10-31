function funval = safeLog(x)

funval = log(1e-100+x);