import os
import re
import variablevals as params
from cc3d import CompuCellSetup
from oral_biofilm3dSteppables import oral_biofilm3dSteppable

# --- AUTO-SYNC XML STEPS ---
xml_file = os.path.join(os.path.dirname(__file__), "oral_biofilm3d.xml")
if os.path.exists(xml_file):
    with open(xml_file, "r") as f:
        xml = f.read()
    xml = re.sub(r"<Steps>\d+</Steps>",
                 f"<Steps>{params.RECOMMENDED_STEPS}</Steps>", xml)
    with open(xml_file, "w") as f:
        f.write(xml)
    print(f"[INFO] Synced XML <Steps> with RECOMMENDED_STEPS={params.RECOMMENDED_STEPS}")
else:
    print(f"[WARNING] XML file {xml_file} not found, skipping sync!")

# --- NORMAL CC3D SETUP ---
try:
    from variablevals import steppable_frequency
except ImportError as err:
    print(f"[Warning] {err}. Using default steppable_frequency=100.")
    steppable_frequency = 100

CompuCellSetup.register_steppable(
    steppable=oral_biofilm3dSteppable(frequency=steppable_frequency)
)
CompuCellSetup.run()