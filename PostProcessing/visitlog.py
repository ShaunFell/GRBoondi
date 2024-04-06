# Visit 3.3.3 log file
ScriptVersion = "3.3.3"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
ShowAllWindows()
OpenDatabase("/home/hd/hd_hd/hd_pb293/WS_GRChombo/testing/GRMilijun_ProcaKerrBH_23287046/hdf5/GeneralizedProcap_*.3d.hdf5 database", 0)
# The UpdateDBPluginInfo RPC is not supported in the VisIt module so it will not be logged.
