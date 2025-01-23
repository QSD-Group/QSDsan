# import qsdsan as qs
# from qsdsan import sanunits as su, processes as pc

# cmps = pc.create_asm1_cmps()
# asm1 = pc.ASM1()
# inf = qs.WasteStream('inf', H2O=1.53e6, S_I=46, S_S=54, X_I=1770, X_S=230, 
#                    X_BH=3870, X_BA=225, X_P=680, S_O=0.377, S_NO=7.98, 
#                    S_NH=25.6, S_ND=5.87, X_ND=13.4, S_ALK=103)
# inf.show()

# eff = qs.WasteStream('eff')
# AS = su.CSTR('O1', ins=inf, outs=eff, DO_ID='S_O', suspended_growth_model=asm1)
# AS.set_init_conc(S_I=30, S_S=5, X_I=1000, X_S=100, X_BH=500, X_BA=100, 
#                       X_P=100, S_O=2, S_NO=20, S_NH=2, S_ND=1, X_ND=1, S_ALK=84)
# AS.simulate(t_span=(0,100))

# eff.show()

import qsdsan.sanunits as su, qsdsan.processes as pc
from qsdsan import WasteStream
cmps = pc.create_asm1_cmps()
asm1 = pc.ASM1()
inf = WasteStream('inf', H2O=1.53e6, S_I=46, S_S=54, X_I=1770, X_S=230, 
                   X_BH=3870, X_BA=225, X_P=680, S_O=0.377, S_NO=7.98, 
                   S_NH=25.6, S_ND=5.87, X_ND=13.4, S_ALK=103)
inf.show()
eff = WasteStream('eff')
AS = su.PFR('AS', ins=inf, outs=eff, tanks_share_walls=True,
             N_tanks_in_series=5, V_tanks=[1000]*2+[1333]*3, 
             influent_fractions=[[1.0, 0,0,0,0]], DO_setpoints=[0]*2+[1.7, 2.4, 0.5],
             internal_recycles=[(4,0,55338)], DO_ID='S_O',
             suspended_growth_model=asm1)
AS.set_init_conc(S_I=30, S_S=5, X_I=1000, X_S=100, X_BH=500, X_BA=100, 
                  X_P=100, S_O=2, S_NO=20, S_NH=2, S_ND=1, X_ND=1, S_ALK=84)
AS.simulate(t_span=(0,100), method='BDF')

eff.show()