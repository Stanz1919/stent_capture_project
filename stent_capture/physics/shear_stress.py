class WallShearStress:
    """
    Class to analyze wall shear stress in vascular systems.
    """

    def __init__(self, fluid_viscosity, flow_rate, pipe_radius):
        """
        Initialize with fluid viscosity, flow rate, and pipe radius.
        :param fluid_viscosity: Viscosity of the fluid (Pa.s)
        :param flow_rate: Volumetric flow rate (m^3/s)
        :param pipe_radius: Radius of the pipe (m)
        """
        self.fluid_viscosity = fluid_viscosity
        self.flow_rate = flow_rate
        self.pipe_radius = pipe_radius

    def calculate_wall_shear_stress(self):
        """
        Calculate the wall shear stress using the formula:
        tau_w = (4 * mu * Q) / (pi * r^3)
        :return: Wall shear stress (Pa)
        """
        import math
        tau_w = (4 * self.fluid_viscosity * self.flow_rate) / (math.pi * (self.pipe_radius ** 3))
        return tau_w


class ShearStressProfile:
    """
    Class to model the shear stress profile along a blood vessel.
    """

    def __init__(self, wall_shear_stress, length):
        """
        Initialize with wall shear stress and length of the vessel.
        :param wall_shear_stress: Wall shear stress (Pa)
        :param length: Length of the vessel (m)
        """
        self.wall_shear_stress = wall_shear_stress
        self.length = length

    def circumferential_profile(self):
        """
        Calculate the circumferential shear stress profile.
        :return: List of circumferential shear stress values
        """
        return [self.wall_shear_stress * (1 - (r/self.length)**2) for r in range(0, int(self.length))]

    def axial_profile(self):
        """
        Calculate the axial shear stress profile.
        :return: List of axial shear stress values
        """
        return [self.wall_shear_stress] * int(self.length)


class CellAttachmentMechanics:
    """
    Class for analyzing cell attachment under shear stress.
    """

    def __init__(self, shear_stress, attachment_strength):
        """
        Initialize with shear stress and cell attachment strength.
        :param shear_stress: Shear stress acting on cells (Pa)
        :param attachment_strength: Strength of cell attachment (Pa)
        """
        self.shear_stress = shear_stress
        self.attachment_strength = attachment_strength

    def shear_force_on_cells(self):
        """
        Calculate the shear force acting on the cells.
        :return: Shear force on cells (N)
        """
        return self.shear_stress * 1.0  # Placeholder for actual surface area

    def detachment_criteria(self):
        """
        Determine if cells will detach based on shear stress and attachment strength.
        :return: Boolean indicating whether cells detach
        """
        return self.shear_stress > self.attachment_strength

    def clinical_parameter_ranges(self):
        """
        Provide clinical parameter ranges for shear stress and attachment strength.
        :return: Dictionary of clinical ranges
        """
        return {
            'shear_stress_range': (0.1, 10.0),  # Example range in Pa
            'attachment_strength_range': (0.5, 5.0)  # Example range in Pa
        }