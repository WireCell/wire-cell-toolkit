# WireCellRoot (mirrors root/wscript_build).  Builds a ROOT dictionary from
# dict/LinkDef.h (wct_package_root_dict hook, wct-ike.7) and the apps/ programs.
wct_package(WireCellRoot
            USE      WireCellAux WireCellClus ROOTSYS
            APP_USE  WireCellAux ROOTSYS
            TEST_USE WireCellGen WireCellSigProc)
