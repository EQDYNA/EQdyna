import netCDF4 as nc
import sys
test = nc.Dataset('test.nc')
ref  = nc.Dataset('ref.nc')

var = sys.argv[1]
diff = test[var][:][:] - ref[var][:][:]

mean_abs = abs(diff).mean()
print('diff mean is', mean_abs)
mean_test = abs(test[var][:][:]).mean()
print('test mean is', mean_test)
