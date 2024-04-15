package WormBot
  package Components
    package Valves
      partial model ValveInterface
        replaceable package Medium = Modelica.Media.Air.DryAirNasa "Medium model";
        Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {-100, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {100, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation(
          Placement(visible = true, transformation(origin = {100, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation(
          Placement(visible = true, transformation(origin = {-100, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        annotation(
          Icon(graphics = {Line(origin = {0.362069, -40.5517}, points = {{-40, 40}, {40, -40}, {40, 40}, {-40, -40}, {-40, 40}}, thickness = 0.5), Rectangle(extent = {{-100, 100}, {100, -100}}), Line(origin = {-69.2379, -40.5517}, points = {{-30, 0}, {30, 0}}, color = {0, 127, 255}), Line(origin = {70.3621, -40.5517}, points = {{-30, 0}, {30, 0}}, color = {0, 127, 255}), Text(origin = {0, -186}, textColor = {0, 0, 255}, extent = {{-150, 85}, {150, 45}}, textString = "%name"), Rectangle(origin = {0, 49}, lineThickness = 0.5, extent = {{-40, 29}, {40, -29}}), Line(origin = {0, -10}, points = {{0, 30}, {0, -30}}, thickness = 0.5), Line(origin = {-70, 60.41}, points = {{-30, 0}, {30, 0}}, color = {0, 0, 243}), Line(origin = {70, 60}, points = {{-30, 0}, {30, 0}}, color = {0, 0, 255})}));
      end ValveInterface;

      model SolenoidValve
        extends WormBot.Components.Valves.ValveInterface;
        parameter Modelica.Units.SI.PressureDifference dp_nominal = 10000 "Nominal pressure difference across the valve";
        parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.0237 "Nominal mass flow rate across the valve at dp_nominal";
        Modelica.Electrical.Analog.Basic.Resistor resistor(R = 12/0.54) annotation(
          Placement(visible = true, transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear2(redeclare package Medium = Medium, dp_nominal(displayUnit = "kPa") = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {0, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression realExpression(y = valve_Open) annotation(
          Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Real valve_Open = (pin_p.i/0.5401);
      equation
        connect(pin_p, resistor.p) annotation(
          Line(points = {{-100, 60}, {-10, 60}}, color = {0, 0, 255}));
        connect(resistor.n, pin_n) annotation(
          Line(points = {{10, 60}, {100, 60}}, color = {0, 0, 255}));
        connect(resistor.p, pin_p) annotation(
          Line(points = {{-10, 60}, {-100, 60}}, color = {0, 0, 255}));
        connect(resistor.n, pin_n) annotation(
          Line(points = {{10, 60}, {100, 60}}, color = {0, 0, 255}));
        connect(valveLinear2.port_a, port_a) annotation(
          Line(points = {{-10, -40}, {-100, -40}}, color = {0, 127, 255}));
        connect(valveLinear2.port_b, port_b) annotation(
          Line(points = {{10, -40}, {100, -40}}, color = {0, 127, 255}));
        connect(realExpression.y, valveLinear2.opening) annotation(
          Line(points = {{-38, 0}, {0, 0}, {0, -32}}, color = {0, 0, 127}));
      end SolenoidValve;
    end Valves;

    package RadialActuators
      partial model RadialInterface
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium model";
        Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a annotation(
          Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b annotation(
          Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        annotation(
          Icon(graphics = {Ellipse(lineColor = {0, 127, 0}, fillColor = {142, 203, 135}, fillPattern = FillPattern.Sphere, extent = {{-60, 60}, {60, -60}}), Line(origin = {-81, 0.27}, points = {{-21, 0}, {21, 0}}, color = {0, 127, 0}), Line(origin = {80, 0}, points = {{20, 0}, {-20, 0}}, color = {0, 127, 0}), Line(origin = {0, 80}, points = {{0, 20}, {0, -20}}), Rectangle(extent = {{-100, 100}, {100, -100}}), Text(origin = {0, -146}, textColor = {0, 0, 255}, extent = {{-150, 85}, {150, 45}}, textString = "%name")}));
      end RadialInterface;

      model RadialPressure
        extends WormBot.Components.RadialActuators.RadialInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Radius r_joint = 0.02 "Joint radius";
        parameter Real gamma_max = 1.5 "Maximum radial stretch ratio";
        parameter Real r_actuator_0 = r_joint*0.9 "Unstretched actuator radius";
        parameter Modelica.Units.SI.Volume V_0 = 0.01 "Unstretched actuator volume";
        parameter Modelica.Units.SI.Mass m = 5;
        parameter Modelica.Units.SI.Pressure p_start = 101325 "Starting pressure";
        parameter Real mu_k_joint = 0.2 "Kinetic friction coefficient for the joint";
        parameter Real mu_s_joint = 0.25 "Static friction coefficient for the joint";
        parameter Real mu_k_actuator = 0.5 "Kinetic friction coefficient for the actuator";
        parameter Real mu_s_actuator = 0.55 "Static friction coefficient for the actuator";
        //
        // Variables
        //
        Modelica.Units.SI.Radius r_actuator "Current radius";
        Real gamma(max = gamma_max, min = 1, start = 1) "Radial stretch ratio";
        Boolean actuatorContact "Is the actuator in contact with the ground?";
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = V_0, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, p_start = p_start, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.ContactForces.OnOffContact onOffContact annotation(
          Placement(visible = true, transformation(origin = {-30, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.BMP frictionJoint(mu_k = mu_k_joint, mu_s = mu_s_joint) annotation(
          Placement(visible = true, transformation(origin = {-20, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.BMP frictionActuator(mu_k = mu_k_actuator, mu_s = mu_s_actuator) annotation(
          Placement(visible = true, transformation(origin = {20, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression normalForce(y = m*9.81) annotation(
          Placement(visible = true, transformation(origin = {-70, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.BooleanExpression actuatorOn(y = actuatorContact) annotation(
          Placement(visible = true, transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
//
// Equations
//
        if r_actuator >= r_joint then
          actuatorContact = true;
        else
          actuatorContact = false;
        end if;
        gamma = pressure.p/101325;
        r_actuator = gamma*r_actuator_0;
//
// Connect equations
//
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {30, 52}, {0, 52}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 52}, {0, 52}, {0, 100}}, color = {0, 127, 255}));
        connect(frictionActuator.flange_a, frictionJoint.flange_b) annotation(
          Line(points = {{10, -30}, {-10, -30}}, color = {0, 127, 0}));
        connect(onOffContact.joint, frictionJoint.fn) annotation(
          Line(points = {{-20, 14}, {0, 14}, {0, -10}, {-20, -10}, {-20, -20}}, color = {0, 0, 127}));
        connect(onOffContact.actuator, frictionActuator.fn) annotation(
          Line(points = {{-20, 26}, {20, 26}, {20, -20}}, color = {0, 0, 127}));
        connect(frictionActuator.flange_b, flange_b) annotation(
          Line(points = {{30, -30}, {80, -30}, {80, 0}, {100, 0}}, color = {0, 127, 0}));
        connect(frictionJoint.flange_a, flange_a) annotation(
          Line(points = {{-30, -30}, {-80, -30}, {-80, 0}, {-100, 0}}, color = {0, 127, 0}));
        connect(normalForce.y, onOffContact.fn) annotation(
          Line(points = {{-59, 40}, {-52, 40}, {-52, 26}, {-40, 26}}, color = {0, 0, 127}));
        connect(actuatorOn.y, onOffContact.switching) annotation(
          Line(points = {{-58, 10}, {-52, 10}, {-52, 14}, {-40, 14}}, color = {255, 0, 255}));
      end RadialPressure;

      model RadialPressureContact
        extends WormBot.Components.RadialActuators.RadialInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Radius r_joint = 0.0385 "Joint radius";
        parameter Modelica.Units.SI.Length l_0 = 0.0625 "Length of radial actuator";
        parameter Real gamma_max = 1.5 "Maximum radial stretch ratio";
        parameter Real r_actuator_0 = 0.0375 "Unstretched actuator radius";
        parameter Modelica.Units.SI.Volume V_0 = Modelica.Constants.pi*l_0*r_actuator_0^(2) "Unstretched actuator volume";
        parameter Modelica.Units.SI.Mass m = 5;
        parameter Modelica.Units.SI.Pressure p_start = 101325 "Starting pressure";
        parameter Real mu_k_joint = 0.2 "Kinetic friction coefficient for the joint";
        parameter Real mu_s_joint = 0.25 "Static friction coefficient for the joint";
        parameter Real mu_k_actuator = 0.5 "Kinetic friction coefficient for the actuator";
        parameter Real mu_s_actuator = 0.55 "Static friction coefficient for the actuator";
        parameter Modelica.Units.SI.Pressure p_ambient = 101325 "Ambient pressure";
        parameter Modelica.Units.SI.ModulusOfElasticity E1 = 58460;//10*10^(6);
        parameter Modelica.Units.SI.ModulusOfElasticity E2 = 10540;        //10*10^(6);
        parameter Modelica.Units.SI.ModulusOfElasticity E = E1+E2 ;                //10*10^(6);
        //
        // Variables
        //
        Modelica.Units.SI.Radius r_actuator(start = r_actuator_0, fixed = true) "Current radius";
        Modelica.Units.SI.Radius a_contact "Contact radius";
        Real gamma(max = gamma_max, min = 1, start = 1) "Radial stretch ratio";
        Boolean actuatorContact "Is the actuator in contact with the ground?";
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = V_0, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, p_start = p_start, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.ContactForces.OnOffContact onOffContact annotation(
          Placement(visible = true, transformation(origin = {-30, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.BMP frictionJoint(mu_k = mu_k_joint, mu_s = mu_s_joint) annotation(
          Placement(visible = true, transformation(origin = {-20, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.BMP frictionActuator(mu_k = mu_k_actuator, mu_s = mu_s_actuator) annotation(
          Placement(visible = true, transformation(origin = {20, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression normalForce(y = m*9.81) annotation(
          Placement(visible = true, transformation(origin = {-70, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.BooleanExpression actuatorOn(y = actuatorContact) annotation(
          Placement(visible = true, transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
//
// Equations
//
        if r_actuator >= r_joint then
          actuatorContact = true;
        else
          actuatorContact = false;
        end if;
        gamma = pressure.p/p_ambient;
        r_actuator = gamma*r_actuator_0;
        frictionActuator.fn = (4/3)*E*a_contact^(3)/r_actuator;
//
// Connect equations
//
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {30, 52}, {0, 52}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 52}, {0, 52}, {0, 100}}, color = {0, 127, 255}));
        connect(frictionActuator.flange_a, frictionJoint.flange_b) annotation(
          Line(points = {{10, -30}, {-10, -30}}, color = {0, 127, 0}));
        connect(onOffContact.joint, frictionJoint.fn) annotation(
          Line(points = {{-20, 14}, {0, 14}, {0, -10}, {-20, -10}, {-20, -20}}, color = {0, 0, 127}));
        connect(onOffContact.actuator, frictionActuator.fn) annotation(
          Line(points = {{-20, 26}, {20, 26}, {20, -20}}, color = {0, 0, 127}));
        connect(frictionActuator.flange_b, flange_b) annotation(
          Line(points = {{30, -30}, {80, -30}, {80, 0}, {100, 0}}, color = {0, 127, 0}));
        connect(frictionJoint.flange_a, flange_a) annotation(
          Line(points = {{-30, -30}, {-80, -30}, {-80, 0}, {-100, 0}}, color = {0, 127, 0}));
        connect(normalForce.y, onOffContact.fn) annotation(
          Line(points = {{-59, 40}, {-52, 40}, {-52, 26}, {-40, 26}}, color = {0, 0, 127}));
        connect(actuatorOn.y, onOffContact.switching) annotation(
          Line(points = {{-58, 10}, {-52, 10}, {-52, 14}, {-40, 14}}, color = {255, 0, 255}));
      end RadialPressureContact;

      model RadialPressureDynamic
        extends WormBot.Components.RadialActuators.RadialInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Radius r_joint = 0.0385 "Joint radius";
        parameter Modelica.Units.SI.Length l_0 = 0.0625 "Length of radial actuator";
        parameter Real r_0 = 0.0375 "Unstretched actuator radius";
        parameter Modelica.Units.SI.Volume V_0 = Modelica.Constants.pi*l_0*r_0^(2) "Unstretched actuator volume";
        parameter Modelica.Units.SI.Mass m = 5;
        parameter Modelica.Units.SI.Pressure p_start = 101325 "Starting pressure";
        parameter Real mu_k_joint = 0.2 "Kinetic friction coefficient for the joint";
        parameter Real mu_s_joint = 0.25 "Static friction coefficient for the joint";
        parameter Real mu_k_actuator = 0.5 "Kinetic friction coefficient for the actuator";
        parameter Real mu_s_actuator = 0.55 "Static friction coefficient for the actuator";
        parameter Modelica.Units.SI.Pressure p_ambient = 101325 "Ambient pressure";
        parameter Modelica.Units.SI.ModulusOfElasticity E1 = 58460;
        parameter Modelica.Units.SI.ModulusOfElasticity E2 = 10540;
        parameter Modelica.Units.SI.Thickness t = 0.005;
        //
        // Variables
        //
        Modelica.Units.SI.Radius r_actuator "Current radius";
        Modelica.Units.SI.Radius a_contact "Contact radius";
        Boolean actuatorContact "Is the actuator in contact with the ground?";
        Real a;
        Real b;
        Real c;
        Modelica.Units.SI.Radius delta_r;
        Modelica.Units.SI.Pressure p "Gauge pressure in the actuator";
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = V_0, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, p_start = p_start, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.ContactForces.OnOffContact onOffContact annotation(
          Placement(visible = true, transformation(origin = {-30, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.BMP frictionJoint(mu_k = mu_k_joint, mu_s = mu_s_joint) annotation(
          Placement(visible = true, transformation(origin = {-20, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.BMP frictionActuator(mu_k = mu_k_actuator, mu_s = mu_s_actuator) annotation(
          Placement(visible = true, transformation(origin = {20, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression normalForce(y = m*9.81) annotation(
          Placement(visible = true, transformation(origin = {-70, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.BooleanExpression actuatorOn(y = actuatorContact) annotation(
          Placement(visible = true, transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
//
// Equations
//
        if r_actuator >= r_joint then
          actuatorContact = true;
        else
          actuatorContact = false;
        end if;
//
// Radius
//
        a = E2*t/(r_0*(r_0 + t));
        b = E1*log((r_0 + t)/r_0) - p;
        c = -r_0*p;
        p = pressure.p - p_ambient;
        delta_r = (-b + sqrt(b^(2) - 4*a*c))/(2*a);
        r_actuator = delta_r + r_0;
        frictionActuator.fn = (4/3)*E1*a_contact^(3)/r_actuator;
//
// Connect equations
//
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {30, 52}, {0, 52}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 52}, {0, 52}, {0, 100}}, color = {0, 127, 255}));
        connect(frictionActuator.flange_a, frictionJoint.flange_b) annotation(
          Line(points = {{10, -30}, {-10, -30}}, color = {0, 127, 0}));
        connect(onOffContact.joint, frictionJoint.fn) annotation(
          Line(points = {{-20, 14}, {0, 14}, {0, -10}, {-20, -10}, {-20, -20}}, color = {0, 0, 127}));
        connect(onOffContact.actuator, frictionActuator.fn) annotation(
          Line(points = {{-20, 26}, {20, 26}, {20, -20}}, color = {0, 0, 127}));
        connect(frictionActuator.flange_b, flange_b) annotation(
          Line(points = {{30, -30}, {80, -30}, {80, 0}, {100, 0}}, color = {0, 127, 0}));
        connect(frictionJoint.flange_a, flange_a) annotation(
          Line(points = {{-30, -30}, {-80, -30}, {-80, 0}, {-100, 0}}, color = {0, 127, 0}));
        connect(normalForce.y, onOffContact.fn) annotation(
          Line(points = {{-59, 40}, {-52, 40}, {-52, 26}, {-40, 26}}, color = {0, 0, 127}));
        connect(actuatorOn.y, onOffContact.switching) annotation(
          Line(points = {{-58, 10}, {-52, 10}, {-52, 14}, {-40, 14}}, color = {255, 0, 255}));
      end RadialPressureDynamic;

      model RadialPressureLS
        extends WormBot.Components.RadialActuators.RadialInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Radius r_joint = 0.0385 "Joint radius";
        parameter Modelica.Units.SI.Length l_0 = 0.0625 "Length of radial actuator";
        parameter Real r_actuator_0 = 0.0375 "Unstretched actuator radius";
        parameter Modelica.Units.SI.Volume V_0 = Modelica.Constants.pi*l_0*r_actuator_0^(2) "Unstretched actuator volume";
        parameter Modelica.Units.SI.Mass m = 5;
        parameter Modelica.Units.SI.Pressure p_start = 101325 "Starting pressure";
        parameter Real mu_k_joint = 0.2 "Kinetic friction coefficient for the joint";
        parameter Real mu_s_joint = 0.25 "Static friction coefficient for the joint";
        parameter Real mu_k_actuator = 0.5 "Kinetic friction coefficient for the actuator";
        parameter Real mu_s_actuator = 0.55 "Static friction coefficient for the actuator";
        parameter Modelica.Units.SI.Pressure p_ambient = 101325 "Ambient pressure";
        parameter Modelica.Units.SI.ModulusOfElasticity E = 10*10^(6);
        parameter Real a = 7.8634*10^(-10);
        parameter Real b = -1.541*10^(-4);
        parameter Real c = 7.6104;
        //
        // Variables
        //
        Modelica.Units.SI.Radius r_actuator(start = r_actuator_0) "Current radius";
        Modelica.Units.SI.Radius a_contact "Contact radius";
        Boolean actuatorContact "Is the actuator in contact with the ground?";
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = V_0, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, p_start = p_start, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.ContactForces.OnOffContact onOffContact annotation(
          Placement(visible = true, transformation(origin = {-30, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.BMP frictionJoint(mu_k = mu_k_joint, mu_s = mu_s_joint) annotation(
          Placement(visible = true, transformation(origin = {-20, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.BMP frictionActuator(mu_k = mu_k_actuator, mu_s = mu_s_actuator) annotation(
          Placement(visible = true, transformation(origin = {20, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression normalForce(y = m*9.81) annotation(
          Placement(visible = true, transformation(origin = {-70, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.BooleanExpression actuatorOn(y = actuatorContact) annotation(
          Placement(visible = true, transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
//
// Equations
//
        if r_actuator >= r_joint then
          actuatorContact = true;
        else
          actuatorContact = false;
        end if;
        r_actuator = a*pressure.p^2 + b*pressure.p + c;
        frictionActuator.fn = (4/3)*E*a_contact^(3)/r_actuator;
//
// Connect equations
//
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {30, 52}, {0, 52}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 52}, {0, 52}, {0, 100}}, color = {0, 127, 255}));
        connect(frictionActuator.flange_a, frictionJoint.flange_b) annotation(
          Line(points = {{10, -30}, {-10, -30}}, color = {0, 127, 0}));
        connect(onOffContact.joint, frictionJoint.fn) annotation(
          Line(points = {{-20, 14}, {0, 14}, {0, -10}, {-20, -10}, {-20, -20}}, color = {0, 0, 127}));
        connect(onOffContact.actuator, frictionActuator.fn) annotation(
          Line(points = {{-20, 26}, {20, 26}, {20, -20}}, color = {0, 0, 127}));
        connect(frictionActuator.flange_b, flange_b) annotation(
          Line(points = {{30, -30}, {80, -30}, {80, 0}, {100, 0}}, color = {0, 127, 0}));
        connect(frictionJoint.flange_a, flange_a) annotation(
          Line(points = {{-30, -30}, {-80, -30}, {-80, 0}, {-100, 0}}, color = {0, 127, 0}));
        connect(normalForce.y, onOffContact.fn) annotation(
          Line(points = {{-59, 40}, {-52, 40}, {-52, 26}, {-40, 26}}, color = {0, 0, 127}));
        connect(actuatorOn.y, onOffContact.switching) annotation(
          Line(points = {{-58, 10}, {-52, 10}, {-52, 14}, {-40, 14}}, color = {255, 0, 255}));
      end RadialPressureLS;
    end RadialActuators;

    package SphericalActuators
      partial model SphereInterface
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium model";
        Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a annotation(
          Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b annotation(
          Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        annotation(
          Icon(graphics = {Ellipse(lineColor = {0, 127, 0}, fillColor = {142, 203, 135}, fillPattern = FillPattern.Sphere, extent = {{-60, 60}, {60, -60}}), Line(origin = {-81, 0.27}, points = {{-21, 0}, {21, 0}}, color = {0, 127, 0}), Line(origin = {80, 0}, points = {{20, 0}, {-20, 0}}, color = {0, 127, 0}), Line(origin = {0, 80}, points = {{0, 20}, {0, -20}}), Rectangle(extent = {{-100, 100}, {100, -100}}), Text(origin = {0, -146}, textColor = {0, 0, 255}, extent = {{-150, 85}, {150, 45}}, textString = "%name")}));
      end SphereInterface;

      model RadialPressure
        extends WormBot.Components.SphericalActuators.SphereInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Radius r_joint = 0.02 "Joint radius";
        parameter Real gamma_max = 1.5 "Maximum radial stretch ratio";
        parameter Real r_actuator_0 = r_joint*0.9 "Unstretched actuator radius";
        parameter Modelica.Units.SI.Volume V_0 = 0.01 "Unstretched actuator volume";
        parameter Modelica.Units.SI.Mass m = 5;
        parameter Modelica.Units.SI.Pressure p_start = 101325 "Starting pressure";
        parameter Real mu_k_joint = 0.2 "Kinetic friction coefficient for the joint";
        parameter Real mu_s_joint = 0.25 "Static friction coefficient for the joint";
        parameter Real mu_k_actuator = 0.5 "Kinetic friction coefficient for the actuator";
        parameter Real mu_s_actuator = 0.55 "Static friction coefficient for the actuator";
        //
        // Variables
        //
        Modelica.Units.SI.Radius r_actuator "Current radius";
        Real gamma(max = gamma_max, min = 1, start = 1) "Radial stretch ratio";
        Boolean actuatorContact "Is the actuator in contact with the ground?";
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = V_0, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, p_start = p_start, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.ContactForces.OnOffContact onOffContact annotation(
          Placement(visible = true, transformation(origin = {-30, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.LuGre frictionJoint(mu_k = mu_k_joint, mu_s = mu_s_joint) annotation(
          Placement(visible = true, transformation(origin = {-20, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.LuGre frictionActuator(mu_k = mu_k_actuator, mu_s = mu_s_actuator) annotation(
          Placement(visible = true, transformation(origin = {20, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression normalForce(y = m*9.81) annotation(
          Placement(visible = true, transformation(origin = {-70, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.BooleanExpression actuatorOn(y = actuatorContact) annotation(
          Placement(visible = true, transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
//
// Equations
//
        if r_actuator >= r_joint then
          actuatorContact = true;
        else
          actuatorContact = false;
        end if;
        gamma = pressure.p/101325;
        r_actuator = gamma*r_actuator_0;
//
// Connect equations
//
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {30, 52}, {0, 52}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 52}, {0, 52}, {0, 100}}, color = {0, 127, 255}));
        connect(frictionActuator.flange_a, frictionJoint.flange_b) annotation(
          Line(points = {{10, -30}, {-10, -30}}, color = {0, 127, 0}));
        connect(onOffContact.joint, frictionJoint.fn) annotation(
          Line(points = {{-20, 14}, {0, 14}, {0, -10}, {-20, -10}, {-20, -20}}, color = {0, 0, 127}));
        connect(onOffContact.actuator, frictionActuator.fn) annotation(
          Line(points = {{-20, 26}, {20, 26}, {20, -20}}, color = {0, 0, 127}));
        connect(frictionActuator.flange_b, flange_b) annotation(
          Line(points = {{30, -30}, {80, -30}, {80, 0}, {100, 0}}, color = {0, 127, 0}));
        connect(frictionJoint.flange_a, flange_a) annotation(
          Line(points = {{-30, -30}, {-80, -30}, {-80, 0}, {-100, 0}}, color = {0, 127, 0}));
        connect(normalForce.y, onOffContact.fn) annotation(
          Line(points = {{-59, 40}, {-52, 40}, {-52, 26}, {-40, 26}}, color = {0, 0, 127}));
        connect(actuatorOn.y, onOffContact.switching) annotation(
          Line(points = {{-58, 10}, {-52, 10}, {-52, 14}, {-40, 14}}, color = {255, 0, 255}));
      end RadialPressure;
    end SphericalActuators;

    package CylindricalActuators
      partial model CylinderInterface
        replaceable package Medium = Modelica.Media.Air.DryAirNasa "Medium model";
        Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a annotation(
          Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b annotation(
          Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        annotation(
          Icon(graphics = {Ellipse(origin = {60, 0}, lineColor = {0, 127, 0}, fillColor = {142, 203, 135}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-20, -40}, {20, 40}}), Rectangle(lineColor = {0, 127, 0}, fillColor = {142, 203, 135}, fillPattern = FillPattern.HorizontalCylinder, lineThickness = 0, extent = {{-60, -40}, {60, 40}}), Ellipse(origin = {-60, 0}, fillColor = {0, 127, 0}, fillPattern = FillPattern.Solid, extent = {{-20, -40}, {20, 40}}), Line(origin = {-90, 0}, points = {{-10, 0}, {10, 0}}, color = {0, 127, 0}), Line(origin = {90, 0}, points = {{10, 0}, {-10, 0}}, color = {0, 127, 0}), Text(origin = {0, -146}, textColor = {0, 0, 255}, extent = {{-150, 85}, {150, 45}}, textString = "%name"), Rectangle(extent = {{-100, 100}, {100, -100}}), Line(origin = {0, 70}, points = {{0, 30}, {0, -30}})}));
      end CylinderInterface;

      model LinearPressure
        extends WormBot.Components.CylindricalActuators.CylinderInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Length L = 0.2 "Unstretched length of the cylinder";
        parameter Real gamma_max = 1.5 "Maximum stretch ratio";
        parameter Modelica.Units.SI.Volume V_0 = 0.01 "Unstretched actuator volume";
        //
        // Variables
        //
        Real gamma(max = gamma_max, min = 1, start = 1) "Stretch ratio";
        Modelica.Units.SI.Length l "Current length";
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = V_0, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        gamma = pressure.p/101325;
        l = gamma*L;
        flange_b.s = flange_a.s + l;
        flange_b.f + flange_a.f = 0;
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {0, 60}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 40}, {0, 40}, {0, 100}}, color = {0, 127, 255}));
        annotation(
          Icon);
      end LinearPressure;

      model LinearActuatorStatic
        extends WormBot.Components.CylindricalActuators.CylinderInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Length L = 0.2 "Unstretched length of the cylinder";
        parameter Modelica.Units.SI.Radius R = 0.1 "Radius of the cylinder (constant)";
        parameter Modelica.Units.SI.Thickness t = 0.001 "Thickness of the actuator material";
        parameter Modelica.Units.SI.ModulusOfElasticity E = 10*10^(6) "Elastic modulus";
        parameter Modelica.Units.SI.Volume V_0 = 0.01 "Unstretched actuator volume";
        parameter Modelica.Units.SI.Pressure p_ref = 101325 "Ambient pressure";
        parameter Real k = 2*L*R^(2)/(E*(2*R*t + t^(2)));
        parameter Modelica.Units.SI.Temperature T_start = 293.15;
        //
        // Variables
        //
        Modelica.Units.SI.Length l(start = L) "Current length";
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, T_start = T_start, V = V_0, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, use_T_start = true, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        l = L + k*(pressure.p - p_ref);
        flange_b.s = flange_a.s + l;
        flange_b.f + flange_a.f = 0;
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {0, 60}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 40}, {0, 40}, {0, 100}}, color = {0, 127, 255}));
        annotation(
          Icon);
      end LinearActuatorStatic;

      model LinearActuatorDynamic
        extends WormBot.Components.CylindricalActuators.CylinderInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Length l_0 = 0.15 "Unstretched length of the cylinder";
        parameter Modelica.Units.SI.Radius r_0 = 0.0375 "Radius of the cylinder (constant)";
        parameter Modelica.Units.SI.Thickness t = 0.001 "Thickness of the actuator material";
        parameter Modelica.Units.SI.ModulusOfElasticity E = 10*10^(6) "Elastic modulus";
        parameter Modelica.Units.SI.Volume V_0 = Modelica.Constants.pi*r_0^2*l_0 "Unstretched actuator volume";
        parameter Modelica.Units.SI.Pressure p_ref = 101325 "Ambient pressure";
        parameter Modelica.Units.SI.Mass m_r = 1 "Mass of radial actuator";
        parameter Modelica.Units.SI.Mass m_l = 1 "Mass of the linear actuator";
        parameter Modelica.Units.SI.ModulusOfElasticity E_1 = 58460 "Modulus";
        parameter Modelica.Units.SI.ModulusOfElasticity E_2 = 10540 "Modulus";
        parameter Modelica.Units.SI.DampingCoefficient b = 1;
        parameter Modelica.Units.SI.Pressure p_ambient = 101325;
        //
        // Variables
        //
        Modelica.Units.SI.Length l(start = l_0, fixed = true) "Current length";
        Modelica.Units.SI.Position x1 "Position of tail actuator";
        Modelica.Units.SI.Position x2 "Position of head actuator";
        Modelica.Units.SI.Velocity v1(start = 0) "Velocity of tail actuator";
        Modelica.Units.SI.Velocity v2(start = 0) "Velocity of head actuator";
        Modelica.Units.SI.Pressure p "Pressure in actuator";
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = V_0, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        flange_b.s = flange_a.s + l;
        flange_a.s = x2;
// Head
        flange_b.s = x1;
// Tail
        v1 = der(x1);
        v2 = der(x2);
        p = pressure.p - p_ambient;
        m_r*der(v1) + (1/4)*(der(v1) + der(v2)) - (V_0/l_0)*(1 + (x2 - x1)/l_0)*(E_1 + E_2*(1 + (x2 - x1)/l_0)) = b*(v2 - v1) + flange_b.f - Modelica.Constants.pi*r_0*p;
        m_r*der(v2) + (1/4)*(der(v1) + der(v2)) + (V_0/l_0)*(1 + (x2 - x1)/l_0)*(E_1 + E_2*(1 + (x2 - x1)/l_0)) = -b*(v2 - v1) + flange_a.f + Modelica.Constants.pi*r_0*p;
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {0, 60}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 40}, {0, 40}, {0, 100}}, color = {0, 127, 255}));
        annotation(
          Icon);
      end LinearActuatorDynamic;

      model LinearActuatorStaticVarVol
        extends WormBot.Components.CylindricalActuators.CylinderInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Length L = 0.2 "Unstretched length of the cylinder";
        parameter Modelica.Units.SI.Radius R = 0.1 "Radius of the cylinder (constant)";
        parameter Modelica.Units.SI.Thickness t = 0.001 "Thickness of the actuator material";
        parameter Modelica.Units.SI.ModulusOfElasticity E = 10*10^(6) "Elastic modulus";
        parameter Modelica.Units.SI.Volume V_0 = 0.01 "Unstretched actuator volume";
        parameter Modelica.Units.SI.Pressure p_ref = 101325 "Ambient pressure";
        parameter Real k = 2*L*R^(2)/(E*(2*R*t + t^(2)));
        parameter Modelica.Units.SI.Temperature T_start = 293.15;
        //
        // Variables
        //
        Modelica.Units.SI.Length l(start = L) "Current length";
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.OpenTank volume(redeclare package Medium = Medium, T_start = T_start, crossArea = 0.1, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, height = 0.1, level(fixed = true, start = 0.01), massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, use_T_start = true, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        l = L + k*(pressure.p - p_ref);
        flange_b.s = flange_a.s + l;
        flange_b.f + flange_a.f = 0;
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {0, 60}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 40}, {0, 40}, {0, 100}}, color = {0, 127, 255}));
        annotation(
          Icon,
          experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
      end LinearActuatorStaticVarVol;

      model LinearActuatorStaticLS
        extends WormBot.Components.CylindricalActuators.CylinderInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Length l_0 = 0.127 "Unstretched length of the cylinder";
        parameter Modelica.Units.SI.Radius r_0 = 0.0375 "Radius of the cylinder (constant)";
        parameter Modelica.Units.SI.Volume V_0 = Modelica.Constants.pi*l_0*r_0^(2) "Unstretched actuator volume";
        parameter Modelica.Units.SI.Pressure p_ref = 101325 "Ambient pressure";
        parameter Real a = 3.6190*10^(-10) "p^2 coefficient"; // -3.5996*10^(-11)
        parameter Real b = 2.1326*10^(-6) "p coefficient"; // 1.5475*10^(-5)
        parameter Real c = 4.4401*10^(-4) "Constant coefficient"; // -1.0461
        parameter Modelica.Units.SI.Temperature T_start = 293.15;
        //
        // Variables
        //
        Modelica.Units.SI.PressureDifference dp(start = 0) "Internal pressure - ambient pressure";
        Modelica.Units.SI.Length l(start = l_0) "Current length";
        Modelica.Units.SI.Length elong;
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, T_start = T_start, V = V_0, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, use_T_start = true, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        dp = pressure.p - p_ref;
//l = a*pressure.p^2 + b*pressure.p + c;
        elong = a*dp^2 + b*dp + c;
        flange_b.s = flange_a.s + l;
        flange_b.f + flange_a.f = 0;
        elong = l - l_0;
//
// Connections
//
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {0, 60}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 40}, {0, 40}, {0, 100}}, color = {0, 127, 255}));
        annotation(
          Icon);
      end LinearActuatorStaticLS;

      model LinearActuatorDynamic1
        extends WormBot.Components.CylindricalActuators.CylinderInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Length l_0 = 0.127 "Unstretched length of the cylinder";
        parameter Modelica.Units.SI.Radius r_0 = 0.0375 "Radius of the cylinder (constant)";
        parameter Modelica.Units.SI.Volume V_0 = Modelica.Constants.pi*l_0*r_0^(2) "Unstretched actuator volume";
        parameter Modelica.Units.SI.Thickness t = 0.001 "Thickness of the actuator material";
        parameter Modelica.Units.SI.Pressure p_ref = 101325 "Ambient pressure";
        parameter Modelica.Units.SI.Mass m_r = 1 "Mass of radial actuator";
        parameter Modelica.Units.SI.Mass m_l = 1 "Mass of the linear actuator";
        parameter Modelica.Units.SI.ModulusOfElasticity E_1 = 58460 "Modulus";
        parameter Modelica.Units.SI.ModulusOfElasticity E_2 = 10540 "Modulus";
        parameter Modelica.Units.SI.DampingCoefficient b = 1;
        parameter Modelica.Units.SI.Temperature T_start = 298;
        //
        // Variables
        //
        Modelica.Units.SI.PressureDifference dp(start = 0) "Internal pressure - ambient pressure";
        Modelica.Units.SI.Length l(start = l_0) "Current length";
        Modelica.Units.SI.Force f1 "Force into the tail actuator";
        Modelica.Units.SI.Force f2 "Force into the head actuator";
        Modelica.Units.SI.Stress sigma "Stress";
        Modelica.Units.SI.Strain epsilon(start = 0) "Strain";
        Modelica.Units.SI.Position x1 "Position of the tail actuator";
        Modelica.Units.SI.Position x2 "Posotion of the head actuator";
        Modelica.Units.SI.Velocity v1(start = 0) "Velocity of tail actuator";
        Modelica.Units.SI.Velocity v2(start = 0) "Velocity of head actuator";
        Modelica.Units.SI.Pressure p "Pressure in the actuator";
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, T_start = T_start, V = V_0, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, use_T_start = true, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
//
// Head
//
        flange_a.s = x2;
        flange_a.f = f2;
        v2 = der(x2);
//
// Tail
//
        flange_b.s = x1;
        flange_b.f = f1;
        v1 = der(x1);
//
// Balance
//
        f1 + f2 = 0;
        l = x1 - x2;
//
//
//
        p = pressure.p;
        epsilon = (l - l_0)/l_0;
        sigma = E_1*epsilon;
        m_r*der(v1) + (1/6)*m_l*(2*der(v1) + der(v2)) - b*(v2 - v1) - (V_0/l_0)*sigma = -Modelica.Constants.pi*r_0^(2)*p + f1;
        m_r*der(v2) + (1/6)*m_l*(2*der(v2) + der(v1)) + b*(v2 - v1) + (V_0/l_0)*sigma = Modelica.Constants.pi*r_0^(2)*p + f2;
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {0, 60}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 40}, {0, 40}, {0, 100}}, color = {0, 127, 255}));
        annotation(
          Icon,
          experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
      end LinearActuatorDynamic1;

      model LinearActuatorDynamicDirection
        extends WormBot.Components.CylindricalActuators.CylinderInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Length l_0 = 0.125 "Unstretched length of the cylinder";
        parameter Modelica.Units.SI.Radius r_0 = 0.0375 "Radius of the cylinder (constant)";
        parameter Modelica.Units.SI.Thickness t = 0.001 "Thickness of the actuator material";
        //parameter Modelica.Units.SI.ModulusOfElasticity E = 10*10^(6) "Elastic modulus";
        parameter Modelica.Units.SI.Volume V_0 = Modelica.Constants.pi*r_0^2*l_0 "Unstretched actuator volume";
        parameter Modelica.Units.SI.Mass m_r = 1 "Mass of radial actuator";
        parameter Modelica.Units.SI.Mass m_l = 1 "Mass of the linear actuator";
        parameter Modelica.Units.SI.ModulusOfElasticity E_1 = 58460 "Modulus";
        parameter Modelica.Units.SI.ModulusOfElasticity E_2 = 10540 "Modulus";
        parameter Modelica.Units.SI.Damping b = 1 "Damping coefficient";
        parameter Modelica.Units.SI.Pressure p_ambient = 101325;
        parameter Real a_poly = -3.5996*10^(-11) "p^2 coefficient";
        parameter Real b_poly = 1.5475*10^(-5) "p coefficient";
        parameter Real c_poly = -1.0461 "Constant coefficient";
        //
        // Variables
        //
        Modelica.Units.SI.Length l(start = l_0, fixed = true) "Current length";
        Modelica.Units.SI.Position x1 "Position of tail actuator";
        Modelica.Units.SI.Position x2 "Position of head actuator";
        Modelica.Units.SI.Velocity v1(start = 0) "Velocity of tail actuator";
        Modelica.Units.SI.Velocity v2(start = 0) "Velocity of head actuator";
        Modelica.Units.SI.Force f1 "Force into the tail actuator";
        Modelica.Units.SI.Force f2 "Force into the head actuator";
        Modelica.Units.SI.Pressure p "Pressure difference in actuator from ambient";
        Modelica.Units.SI.Force a1;
        Modelica.Units.SI.Force a2;
        Modelica.Units.SI.Force a3;
        Modelica.Units.SI.Force a4;
        Modelica.Units.SI.Force a5;
        Modelica.Units.SI.Force a6;
        Modelica.Units.SI.Force a7;
        Real strain(start = 0);
        Modelica.Units.SI.Stress stress(start = 0);
        Modelica.Units.SI.Length elong;
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = V_0, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        l = x1 - x2;
        flange_a.s = x2;
        flange_a.f = f2;
// Head
        flange_b.s = x1;
        flange_b.f = f1;
// Tail
        v1 = der(x1);
        v2 = der(x2);
        p = pressure.p - p_ambient;
//strain = (a_poly*pressure.p^2 + b_poly*pressure.p + c_poly - l_0)/l_0;
        strain = (l - l_0)/l_0;
        stress = (E_1+E_2)*strain; //strain*(E_1 + E_2*strain);
        a1 = m_r*der(v1);
        a2 = (1/6)*m_l*(2*der(v1) + der(v2));
        a3 = (V_0/l_0)*stress;
        a4 = b*(der(l));
        a5 = Modelica.Constants.pi*r_0^(2)*p;
        a6 = m_r*der(v2);
        a7 = (1/6)*m_l*(2*der(v2) + der(v1));
//
        a1 + a2 + a3 = -a4 + f1 + a5;
        a6 + a7 - a3 = a4 + f2 - a5;
//
        elong = l - l_0;
//
// Connections
//
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {0, 60}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 40}, {0, 40}, {0, 100}}, color = {0, 127, 255}));
        annotation(
          Icon,
          experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
      end LinearActuatorDynamicDirection;
      
      model LinearActuatorDynamicStrain
        extends WormBot.Components.CylindricalActuators.CylinderInterface;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Length l_0 = 0.125 "Unstretched length of the cylinder";
        parameter Modelica.Units.SI.Radius r_0 = 0.0375 "Radius of the cylinder (constant)";
        parameter Modelica.Units.SI.Thickness t = 0.001 "Thickness of the actuator material";
        //parameter Modelica.Units.SI.ModulusOfElasticity E = 10*10^(6) "Elastic modulus";
        parameter Modelica.Units.SI.Volume V_0 = Modelica.Constants.pi*r_0^2*l_0 "Unstretched actuator volume";
        parameter Modelica.Units.SI.Mass m_r = 1 "Mass of radial actuator";
        parameter Modelica.Units.SI.Mass m_l = 1 "Mass of the linear actuator";
        parameter Modelica.Units.SI.ModulusOfElasticity E_1 = 58460 "Modulus";
        parameter Modelica.Units.SI.ModulusOfElasticity E_2 = 10540 "Modulus";
        parameter Modelica.Units.SI.Damping b = 1 "Damping coefficient";
        parameter Modelica.Units.SI.Pressure p_ambient = 101325;
        parameter Real a_poly = 3.6190*10^(-10) "p^2 coefficient"; // -3.5996*10^(-11)
        parameter Real b_poly = 2.1326*10^(-6) "p coefficient"; // 1.5475*10^(-5)
        parameter Real c_poly = 4.4401*10^(-4) "Constant coefficient";         // -1.0461
        //parameter Real a_poly = -3.5996*10^(-11) "p^2 coefficient";
        //parameter Real b_poly = 1.5475*10^(-5) "p coefficient";
        //parameter Real c_poly = -1.0461 "Constant coefficient";
        //
        // Variables
        //
        Modelica.Units.SI.Length l(start = l_0, fixed = true) "Current length";
        Modelica.Units.SI.Position x1 "Position of tail actuator";
        Modelica.Units.SI.Position x2 "Position of head actuator";
        Modelica.Units.SI.Velocity v1(start = 0) "Velocity of tail actuator";
        Modelica.Units.SI.Velocity v2(start = 0) "Velocity of head actuator";
        Modelica.Units.SI.Force f1 "Force into the tail actuator";
        Modelica.Units.SI.Force f2 "Force into the head actuator";
        Modelica.Units.SI.Pressure p "Pressure difference in actuator from ambient";
        Modelica.Units.SI.Force a1;
        Modelica.Units.SI.Force a2;
        Modelica.Units.SI.Force a3;
        Modelica.Units.SI.Force a4;
        Modelica.Units.SI.Force a5;
        Modelica.Units.SI.Force a6;
        Modelica.Units.SI.Force a7;
        Real strain(start = 0);
        Modelica.Units.SI.Stress stress(start = 0);
        Modelica.Units.SI.Length elong;
        //
        // Components
        //
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = V_0, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, nPorts = 1, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        l = x1 - x2;
        flange_a.s = x2;
        flange_a.f = f2;
// Head
        flange_b.s = x1;
        flange_b.f = f1;
// Tail
        v1 = der(x1);
        v2 = der(x2);
        p = pressure.p - p_ambient;
//strain = (a_poly*pressure.p^2 + b_poly*pressure.p + c_poly - l_0)/l_0;
        strain = (elong)/l_0;
        elong = a_poly*p^2 + b_poly*p + c_poly;
        stress = (62150)*strain^(0.7036); //strain*(E_1 + E_2*strain);
        a1 = m_r*der(v1);
        a2 = (1/6)*m_l*(2*der(v1) + der(v2));
        a3 = (V_0/l_0)*stress;
        a4 = b*(der(l));
        a5 = Modelica.Constants.pi*r_0^(2)*p;
        a6 = m_r*der(v2);
        a7 = (1/6)*m_l*(2*der(v2) + der(v1));
//
        a1 + a2 + a3 = -a4 + f1 + a5;
        a6 + a7 - a3 = a4 + f2 - a5;
//
        elong = l - l_0;
//
// Connections
//
        connect(pressure.port, port_a) annotation(
          Line(points = {{30, 60}, {0, 60}, {0, 100}}, color = {0, 127, 255}));
        connect(volume.ports[1], port_a) annotation(
          Line(points = {{-48, 60}, {-48, 40}, {0, 40}, {0, 100}}, color = {0, 127, 255}));
        annotation(
          Icon,
          experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
      end LinearActuatorDynamicStrain;
    end CylindricalActuators;

    package ValveControl
      partial model ControlInterface
        Modelica.Electrical.Analog.Interfaces.PositivePin pin_p[6] annotation(
          Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin_n[6] annotation(
          Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        annotation(
          Diagram,
          Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Rectangle(origin = {0, 20}, lineColor = {95, 95, 95}, fillColor = {156, 156, 156}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-60, 40}, {60, -40}}), Rectangle(origin = {0, -30}, lineColor = {95, 95, 95}, fillColor = {156, 156, 156}, fillPattern = FillPattern.VerticalCylinder, extent = {{-10, -10}, {10, 10}}), Rectangle(origin = {0, -44}, lineColor = {95, 95, 95}, fillColor = {95, 95, 95}, fillPattern = FillPattern.Solid, extent = {{-60, 4}, {60, -4}}), Rectangle(origin = {0, 20}, lineColor = {95, 95, 95}, fillColor = {85, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-58, 38}, {58, -38}}), Text(origin = {0, -186}, textColor = {0, 0, 255}, extent = {{-150, 85}, {150, 45}}, textString = "%name")}));
      end ControlInterface;

      model LUTControl
        extends WormBot.Components.ValveControl.ControlInterface;
        Modelica.Blocks.Sources.ContinuousClock clock(startTime = 0) annotation(
          Placement(visible = true, transformation(origin = {-70, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage[6](each V = 12) annotation(
          Placement(visible = true, transformation(origin = {-30, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Tables.CombiTable1Ds controllerLUT(columns = 2:7, extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint, fileName = "C:/Users/adams/OneDrive - Carleton University/Desktop/Adam Files/Winter 2024/Systems Modelling/Project/Repo/WormBot/WormControl.txt", smoothness = Modelica.Blocks.Types.Smoothness.ConstantSegments, tableName = "wormControl", tableOnFile = true, verboseRead = false) annotation(
          Placement(visible = true, transformation(origin = {-30, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Electrical.Analog.Ideal.IdealClosingSwitch switch[6](each Goff = 1E-9, each Ron = 1E-9) annotation(
          Placement(visible = true, transformation(origin = {30, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Logical.GreaterThreshold greaterThreshold[6](each threshold = 0.1) annotation(
          Placement(visible = true, transformation(origin = {10, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
// Get the control from the LUT
//
        connect(clock.y, controllerLUT.u) annotation(
          Line(points = {{-58, 80}, {-42, 80}}, color = {0, 0, 127}));
        connect(pin_n, constantVoltage.n) annotation(
          Line(points = {{-100, 0}, {-40, 0}}, color = {0, 0, 255}, thickness = 0.5));
        connect(constantVoltage.p, switch.p) annotation(
          Line(points = {{-20, 0}, {20, 0}}, color = {0, 0, 255}, thickness = 0.5));
        connect(switch.n, pin_p) annotation(
          Line(points = {{40, 0}, {100, 0}}, color = {0, 0, 255}, thickness = 0.5));
        connect(controllerLUT.y, greaterThreshold.u) annotation(
          Line(points = {{-18, 80}, {-2, 80}}, color = {0, 0, 127}, thickness = 0.5));
        connect(greaterThreshold.y, switch.control) annotation(
          Line(points = {{22, 80}, {30, 80}, {30, 12}}, color = {255, 0, 255}, thickness = 0.5));
      end LUTControl;
    end ValveControl;

    package ContactForces
      model OnOffContact
        //
        // Inputs
        //
        //
        // Parameters
        //
        //
        // Variables
        //
        Modelica.Blocks.Logical.Switch switchActuator annotation(
          Placement(visible = true, transformation(origin = {30, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Logical.Switch switchJoint annotation(
          Placement(visible = true, transformation(origin = {30, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression nearZero(y = 10^(-10)) annotation(
          Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Interfaces.RealInput fn annotation(
          Placement(visible = true, transformation(origin = {-100, 62}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 62}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        Modelica.Blocks.Interfaces.BooleanInput switching annotation(
          Placement(visible = true, transformation(origin = {-100, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        Modelica.Blocks.Interfaces.RealOutput actuator annotation(
          Placement(visible = true, transformation(origin = {100, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Interfaces.RealOutput joint annotation(
          Placement(visible = true, transformation(origin = {100, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(nearZero.y, switchActuator.u3) annotation(
          Line(points = {{-9, 0}, {0, 0}, {0, 6}, {18, 6}}, color = {0, 0, 127}));
        connect(switchActuator.u1, fn) annotation(
          Line(points = {{18, 22}, {-60, 22}, {-60, 62}, {-100, 62}}, color = {0, 0, 127}));
        connect(switchJoint.u3, fn) annotation(
          Line(points = {{18, -22}, {-60, -22}, {-60, 62}, {-100, 62}}, color = {0, 0, 127}));
        connect(switching, switchJoint.u2) annotation(
          Line(points = {{-100, -60}, {-40, -60}, {-40, -14}, {18, -14}}, color = {255, 0, 255}));
        connect(switching, switchActuator.u2) annotation(
          Line(points = {{-100, -60}, {-40, -60}, {-40, 14}, {18, 14}}, color = {255, 0, 255}));
        connect(switchActuator.y, actuator) annotation(
          Line(points = {{42, 14}, {80, 14}, {80, 60}, {100, 60}}, color = {0, 0, 127}));
        connect(switchJoint.y, joint) annotation(
          Line(points = {{42, -14}, {80, -14}, {80, -60}, {100, -60}}, color = {0, 0, 127}));
        connect(nearZero.y, switchJoint.u1) annotation(
          Line(points = {{-8, 0}, {0, 0}, {0, -6}, {18, -6}}, color = {0, 0, 127}));
      end OnOffContact;
    end ContactForces;

    package Friction
      model LuGre
        //
        // Parameters
        //
        parameter Real mu_k = 0.5;
        parameter Real mu_s = 0.55;
        parameter Real sigma0 = 946160;
        parameter Real sigma1 = 0;
        parameter Real sigma2 = 0;
        parameter Real deltavs = 0.791;
        parameter Real vs = 1.646*10^(-3);
        //
        // Variables
        //
        Modelica.Units.SI.Force f "Friction force";
        Modelica.Units.SI.Force fc "Coulomb force";
        Modelica.Units.SI.Force fs "Static force";
        Modelica.Units.SI.Position s "position";
        Modelica.Units.SI.Velocity v(start = 0, fixed = true) "Sliding velocity";
        Modelica.Units.SI.Acceleration a "Sliding acceleration";
        Real z(start = 0, fixed = true) "Bristle deflection";
        //
        // Components
        //
        Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a annotation(
          Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b annotation(
          Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Interfaces.RealInput fn annotation(
          Placement(visible = true, transformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
      protected
        Modelica.Units.SI.Position s_support = 0 "Absolute postition of the support";
      equation
        flange_a.f + flange_b.f - f = 0;
// Friction
        fc = mu_k*fn;
        fs = mu_s*fn;
        s = flange_a.s - s_support;
        flange_a.s = flange_b.s;
        v = der(s);
        a = der(v);
//
        der(z) = v - (sigma0*abs(v)/(fc + (fs - fc)*exp(-abs(v/vs)^(deltavs))))*z;
        f = sigma0*z + sigma1*der(z) + sigma2*v;
        annotation(
          Icon(graphics = {Text(origin = {0, -146}, textColor = {0, 0, 255}, extent = {{-150, 85}, {150, 45}}, textString = "%name"), Rectangle(extent = {{-100, 100}, {100, -100}}), Line(points = {{-100, 0}, {100, 0}}, color = {0, 127, 0}), Rectangle(lineColor = {0, 127, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-40, 20}, {40, -20}}), Line(origin = {0, -20}, points = {{-60, 0}, {60, 0}}), Rectangle(origin = {0, -26}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Backward, extent = {{-60, 6}, {60, -6}})}));
      end LuGre;

      model Sticking
        //
        // Parameters
        //
        parameter Real mu_k = 0.5 "Coefficient of kinemetic friction.";
        parameter Real mu_s = 0.55 "Coefficient of static friction.";
        //
        // Variables
        //
        Modelica.Units.SI.Force f "Friction force";
        Modelica.Units.SI.Force fc "Coulomb force";
        Modelica.Units.SI.Force fs "Static force";
        Modelica.Units.SI.Force fa "External force at a";
        Modelica.Units.SI.Position s "position";
        Modelica.Units.SI.Velocity v(start = 0, fixed = true) "Sliding velocity";
        Modelica.Units.SI.Acceleration a "Sliding acceleration";
        //
        // Components
        //
        Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a annotation(
          Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b annotation(
          Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Interfaces.RealInput fn annotation(
          Placement(visible = true, transformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
      protected
        Modelica.Units.SI.Position s_support = 0 "Absolute postition of the support";
      equation
        flange_a.f + flange_b.f - f = 0;
// Friction
        fa = flange_a.f + flange_b.f;
        fc = mu_k*fn;
        fs = mu_s*fn;
        s = flange_a.s - s_support;
        flange_a.s = flange_b.s;
        v = der(s);
        a = der(v);
//
        if abs(v) > 1E-7 and abs(fa) < fs then
          f = fa;
        else
          f = sign(v)*fc;
        end if;
        annotation(
          Icon(graphics = {Text(origin = {0, -146}, textColor = {0, 0, 255}, extent = {{-150, 85}, {150, 45}}, textString = "%name"), Rectangle(extent = {{-100, 100}, {100, -100}}), Line(points = {{-100, 0}, {100, 0}}, color = {0, 127, 0}), Rectangle(lineColor = {0, 127, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-40, 20}, {40, -20}}), Line(origin = {0, -20}, points = {{-60, 0}, {60, 0}}), Rectangle(origin = {0, -26}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Backward, extent = {{-60, 6}, {60, -6}})}));
      end Sticking;

      model BMP
        //
        // Parameters
        //
        parameter Real mu_k = 0.5;
        parameter Real mu_s = 0.55;
        parameter Modelica.Units.SI.Velocity v_t = 0.0001 "Transition velocity"; //0.0001
        parameter Modelica.Units.SI.Force fnt = 1 "Transition force";
        //
        // Variables
        //
        Modelica.Units.SI.Force f "Friction force";
        Modelica.Units.SI.Force fd "Dynamic force";
        Modelica.Units.SI.Force fs "Static force";
        Modelica.Units.SI.Position s "position";
        Modelica.Units.SI.Velocity v(start = 0, fixed = true) "Sliding velocity";
        Modelica.Units.SI.Acceleration a "Sliding acceleration";
        //
        // Components
        //
        Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a annotation(
          Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b annotation(
          Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Interfaces.RealInput fn annotation(
          Placement(visible = true, transformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
      protected
        Modelica.Units.SI.Position s_support = 0 "Absolute postition of the support";
      equation
        flange_a.f + flange_b.f - f = 0;
// Friction
        s = flange_a.s - s_support;
        flange_a.s = flange_b.s;
        v = der(s);
        a = der(v);
//
        fd = fn*mu_k*tanh(4*v/v_t);
        fs = fn*(mu_s - mu_k)*(v/v_t)/(0.25*(v/v_t)^(2) + 0.75)^2;
        f = fd + fs;
        annotation(
          Icon(graphics = {Text(origin = {0, -146}, textColor = {0, 0, 255}, extent = {{-150, 85}, {150, 45}}, textString = "%name"), Rectangle(extent = {{-100, 100}, {100, -100}}), Line(points = {{-100, 0}, {100, 0}}, color = {0, 127, 0}), Rectangle(lineColor = {0, 127, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-40, 20}, {40, -20}}), Line(origin = {0, -20}, points = {{-60, 0}, {60, 0}}), Rectangle(origin = {0, -26}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Backward, extent = {{-60, 6}, {60, -6}})}));
      end BMP;
    end Friction;

    model VariableVolume "Volume of fixed size, closed to the ambient, with inlet/outlet ports"
      import Modelica.Constants.pi;
      // Mass and energy balance, ports
      extends Modelica.Fluid.Vessels.BaseClasses.PartialLumpedVessel(final fluidVolume = V, vesselArea = pi*(3/4*V)^(2/3), heatTransfer(surfaceAreas = {4*pi*(3/4*V/pi)^(2/3)}));
      Modelica.Units.SI.Volume V "Volume";
    equation
      Wb_flow = 0;
      for i in 1:nPorts loop
        vessel_ps_static[i] = medium.p;
      end for;
      V = 0.01 + 1E-7*ports[1].p;
      annotation(
        defaultComponentName = "volume",
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Ellipse(extent = {{-100, 100}, {100, -100}}, fillPattern = FillPattern.Sphere, fillColor = {170, 213, 255}), Text(extent = {{-150, 12}, {150, -18}}, textString = "V=%V")}),
        Documentation(info = "<html>
    <p>
    Ideally mixed volume of constant size with two fluid ports and one medium model.
    The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model if <code>use_portsData=false</code>.
    Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected.
    A spherical shape is assumed for the heat transfer area, with V=4/3*pi*r^3, A=4*pi*r^2.
    Ideal heat transfer is assumed per default; the thermal port temperature is equal to the medium temperature.
    </p>
    <p>
    If <code>use_portsData=true</code>, the port pressures represent the pressures just after the outlet (or just before the inlet) in the attached pipe.
    The hydraulic resistances <code>portsData.zeta_in</code> and <code>portsData.zeta_out</code> determine the dissipative pressure drop between volume and port depending on
    the direction of mass flow. See <a href=\"modelica://Modelica.Fluid.Vessels.BaseClasses.VesselPortsData\">VesselPortsData</a> and <em>[Idelchik, Handbook of Hydraulic Resistance, 2004]</em>.
    </p>
    </html>"));
    end VariableVolume;

    model AirIn
      replaceable package Medium = Modelica.Media.Air.DryAirNasa;
      parameter Modelica.Units.SI.MassFlowRate m_flow_regulator;
      input Real onOff[3];
      Modelica.Blocks.Math.Gain gain[3](each k = m_flow_regulator) annotation(
        Placement(visible = true, transformation(origin = {-24, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.MassFlowSource_T pressurizedAir[3](redeclare package Medium = Medium, each nPorts = 1, each use_m_flow_in = true) annotation(
        Placement(visible = true, transformation(origin = {30, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b port_b[3](redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      gain[1].u = onOff[1];
      gain[2].u = onOff[2];
      gain[3].u = onOff[3];
      connect(gain.y, pressurizedAir.m_flow_in) annotation(
        Line(points = {{-13, 0}, {-5, 0}, {-5, 8}, {19, 8}}, color = {0, 0, 127}, thickness = 0.5));
      connect(pressurizedAir[:].ports[1], port_b[:]) annotation(
        Line(points = {{40, 0}, {100, 0}}, color = {0, 127, 255}, thickness = 0.5));
      annotation(
        Icon(graphics = {Text(textColor = {0, 0, 255}, extent = {{-150, 110}, {150, 150}}, textString = "%name"), Ellipse(fillColor = {0, 127, 255}, fillPattern = FillPattern.Sphere, extent = {{-100, 100}, {100, -100}})}));
    end AirIn;
  end Components;

  package Systems
    partial model WormTemplate
      replaceable package Medium = Modelica.Media.Air.DryAirNasa "Medium model";
      WormBot.Components.SphericalActuators.SphereInterface head annotation(
        Placement(visible = true, transformation(origin = {-100, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.SphericalActuators.SphereInterface tail annotation(
        Placement(visible = true, transformation(origin = {100, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      replaceable WormBot.Components.CylindricalActuators.CylinderInterface body annotation(
        Placement(visible = true, transformation(origin = {0, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      replaceable WormBot.Components.Valves.ValveInterface valveHeadIn annotation(
        Placement(visible = true, transformation(origin = {-120, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable WormBot.Components.Valves.ValveInterface valveHeadOut annotation(
        Placement(visible = true, transformation(origin = {-80, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable WormBot.Components.Valves.ValveInterface valveBodyOut annotation(
        Placement(visible = true, transformation(origin = {20, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable WormBot.Components.Valves.ValveInterface valveBodyIn annotation(
        Placement(visible = true, transformation(origin = {-20, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable WormBot.Components.Valves.ValveInterface valveTailOut annotation(
        Placement(visible = true, transformation(origin = {120, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable WormBot.Components.Valves.ValveInterface valveTailIn annotation(
        Placement(visible = true, transformation(origin = {80, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable WormBot.Components.ValveControl.ControlInterface control annotation(
        Placement(visible = true, transformation(origin = {-130, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3) annotation(
        Placement(visible = true, transformation(origin = {160, 100}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 3, p = Medium.p_default*2) annotation(
        Placement(visible = true, transformation(origin = {-180, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-170, 130}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
        Placement(visible = true, transformation(origin = {-170, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(tail.flange_a, body.flange_b) annotation(
        Line(points = {{80, -40}, {20, -40}}, color = {0, 127, 0}));
      connect(body.flange_a, head.flange_b) annotation(
        Line(points = {{-20, -40}, {-80, -40}}, color = {0, 127, 0}));
      connect(valveHeadIn.port_b, head.port_a) annotation(
        Line(points = {{-110, 50}, {-100, 50}, {-100, -20}}, color = {0, 127, 255}));
      connect(head.port_a, valveHeadOut.port_a) annotation(
        Line(points = {{-100, -20}, {-100, 18}, {-90, 18}}, color = {0, 127, 255}));
      connect(valveBodyIn.port_b, body.port_a) annotation(
        Line(points = {{-10, 50}, {0, 50}, {0, -20}}, color = {0, 127, 255}));
      connect(body.port_a, valveBodyOut.port_a) annotation(
        Line(points = {{0, -20}, {0, 20}, {10, 20}}, color = {0, 127, 255}));
      connect(valveTailIn.port_b, tail.port_a) annotation(
        Line(points = {{90, 50}, {100, 50}, {100, -20}}, color = {0, 127, 255}));
      connect(tail.port_a, valveTailOut.port_a) annotation(
        Line(points = {{100, -20}, {100, 20}, {110, 20}}, color = {0, 127, 255}));
      connect(pressurizedAir.ports[1], valveHeadIn.port_a) annotation(
        Line(points = {{-170, 80}, {-140, 80}, {-140, 50}, {-130, 50}}, color = {0, 127, 255}));
      connect(pressurizedAir.ports[2], valveBodyIn.port_a) annotation(
        Line(points = {{-170, 80}, {-40, 80}, {-40, 50}, {-30, 50}}, color = {0, 127, 255}));
      connect(pressurizedAir.ports[3], valveTailIn.port_a) annotation(
        Line(points = {{-170, 80}, {60, 80}, {60, 50}, {70, 50}}, color = {0, 127, 255}));
      connect(valveHeadOut.port_b, ambient.ports[1]) annotation(
        Line(points = {{-70, 18}, {-60, 18}, {-60, 100}, {150, 100}}, color = {0, 127, 255}));
      connect(valveBodyOut.port_b, ambient.ports[2]) annotation(
        Line(points = {{30, 20}, {40, 20}, {40, 100}, {150, 100}}, color = {0, 127, 255}));
      connect(valveTailOut.port_b, ambient.ports[3]) annotation(
        Line(points = {{130, 20}, {140, 20}, {140, 100}, {150, 100}}, color = {0, 127, 255}));
      connect(control.pin_n, ground.p) annotation(
        Line(points = {{-140, 130}, {-160, 130}}, color = {0, 0, 255}, thickness = 0.5));
      connect(control.pin_p[1], valveHeadIn.pin_p) annotation(
        Line(points = {{-120, 130}, {-110, 130}, {-110, 100}, {-136, 100}, {-136, 56}, {-130, 56}}, color = {0, 0, 255}));
      connect(control.pin_p[2], valveHeadOut.pin_p) annotation(
        Line(points = {{-120, 130}, {-96, 130}, {-96, 24}, {-90, 24}}, color = {0, 0, 255}));
      connect(control.pin_p[3], valveBodyIn.pin_p) annotation(
        Line(points = {{-120, 130}, {-36, 130}, {-36, 56}, {-30, 56}}, color = {0, 0, 255}));
      connect(control.pin_p[4], valveBodyOut.pin_p) annotation(
        Line(points = {{-120, 130}, {4, 130}, {4, 26}, {10, 26}}, color = {0, 0, 255}));
      connect(control.pin_p[5], valveTailIn.pin_p) annotation(
        Line(points = {{-120, 130}, {64, 130}, {64, 56}, {70, 56}}, color = {0, 0, 255}));
      connect(control.pin_p[6], valveTailOut.pin_p) annotation(
        Line(points = {{-120, 130}, {104, 130}, {104, 26}, {110, 26}}, color = {0, 0, 255}));
      connect(valveHeadIn.pin_n, ground[1].p) annotation(
        Line(points = {{-110, 56}, {-100, 56}, {-100, 110}, {-160, 110}, {-160, 130}}, color = {0, 0, 255}));
      connect(valveHeadOut.pin_n, ground[2].p) annotation(
        Line(points = {{-70, 24}, {-64, 24}, {-64, 110}, {-160, 110}, {-160, 130}}, color = {0, 0, 255}));
      connect(valveBodyIn.pin_n, ground[3].p) annotation(
        Line(points = {{-10, 56}, {0, 56}, {0, 110}, {-160, 110}, {-160, 130}}, color = {0, 0, 255}));
      connect(valveBodyOut.pin_n, ground[4].p) annotation(
        Line(points = {{30, 26}, {36, 26}, {36, 110}, {-160, 110}, {-160, 130}}, color = {0, 0, 255}));
      connect(valveTailIn.pin_n, ground[5].p) annotation(
        Line(points = {{90, 56}, {100, 56}, {100, 110}, {-160, 110}, {-160, 130}}, color = {0, 0, 255}));
      connect(valveTailOut.pin_n, ground[6].p) annotation(
        Line(points = {{130, 26}, {136, 26}, {136, 110}, {-160, 110}, {-160, 130}}, color = {0, 0, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-200, 140}, {180, -60}})));
    end WormTemplate;

    model WormSystem1
      extends WormBot.Systems.WormTemplate(redeclare replaceable WormBot.Components.ValveControl.LUTControl control, redeclare replaceable WormBot.Components.Valves.SolenoidValve valveHeadIn(redeclare package Medium = Medium), redeclare replaceable WormBot.Components.Valves.SolenoidValve valveHeadOut(redeclare package Medium = Medium), redeclare replaceable WormBot.Components.Valves.SolenoidValve valveBodyIn(redeclare package Medium = Medium), redeclare replaceable WormBot.Components.Valves.SolenoidValve valveBodyOut(redeclare package Medium = Medium), redeclare replaceable WormBot.Components.Valves.SolenoidValve valveTailIn(redeclare package Medium = Medium), redeclare replaceable WormBot.Components.Valves.SolenoidValve valveTailOut(redeclare package Medium = Medium), redeclare replaceable WormBot.Components.SphericalActuators.RadialPressure head(redeclare package Medium = Medium), redeclare replaceable WormBot.Components.CylindricalActuators.LinearActuatorDynamic body(redeclare package Medium = Medium), redeclare replaceable WormBot.Components.SphericalActuators.RadialPressure tail(redeclare package Medium = Medium), system.energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, system.massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial, system.momentumDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial);
      Modelica.Mechanics.Translational.Components.Mass mass(m = 1) annotation(
        Placement(visible = true, transformation(origin = {-170, -40}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(mass.flange_a, head.flange_a) annotation(
        Line(points = {{-160, -40}, {-120, -40}}, color = {0, 127, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 671, Tolerance = 1e-06, Interval = 1.3501));
    end WormSystem1;

    model WormSystem2
      replaceable package Medium = Modelica.Media.Air.DryAirNasa "Medium model";
      WormBot.Components.SphericalActuators.RadialPressure head(redeclare package Medium = Medium, gamma_max = 2, p_start = system.p_start) annotation(
        Placement(visible = true, transformation(origin = {-100, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.SphericalActuators.RadialPressure tail(redeclare package Medium = Medium, gamma_max = 2, p_start = system.p_start) annotation(
        Placement(visible = true, transformation(origin = {100, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      replaceable WormBot.Components.CylindricalActuators.LinearActuatorStatic body(redeclare package Medium = Medium, L = 0.150, R = 0.0375, t = 0.002) annotation(
        Placement(visible = true, transformation(origin = {0, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      replaceable WormBot.Components.Valves.SolenoidValve valveHeadIn(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-120, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable WormBot.Components.Valves.SolenoidValve valveHeadOut(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-80, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable WormBot.Components.Valves.SolenoidValve valveBodyOut(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {20, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable WormBot.Components.Valves.SolenoidValve valveBodyIn(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-20, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable WormBot.Components.Valves.SolenoidValve valveTailOut(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {120, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable WormBot.Components.Valves.SolenoidValve valveTailIn(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {80, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable WormBot.Components.ValveControl.LUTControl control annotation(
        Placement(visible = true, transformation(origin = {-130, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3) annotation(
        Placement(visible = true, transformation(origin = {160, 100}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 3, p = Medium.p_default*1.2) annotation(
        Placement(visible = true, transformation(origin = {-180, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-170, 130}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, momentumDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
        Placement(visible = true, transformation(origin = {-170, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(tail.flange_a, body.flange_b) annotation(
        Line(points = {{80, -40}, {20, -40}}, color = {0, 127, 0}));
      connect(body.flange_a, head.flange_b) annotation(
        Line(points = {{-20, -40}, {-80, -40}}, color = {0, 127, 0}));
      connect(valveHeadIn.port_b, head.port_a) annotation(
        Line(points = {{-110, 50}, {-100, 50}, {-100, -20}}, color = {0, 127, 255}));
      connect(head.port_a, valveHeadOut.port_a) annotation(
        Line(points = {{-100, -20}, {-100, 18}, {-90, 18}}, color = {0, 127, 255}));
      connect(valveBodyIn.port_b, body.port_a) annotation(
        Line(points = {{-10, 50}, {0, 50}, {0, -20}}, color = {0, 127, 255}));
      connect(body.port_a, valveBodyOut.port_a) annotation(
        Line(points = {{0, -20}, {0, 20}, {10, 20}}, color = {0, 127, 255}));
      connect(valveTailIn.port_b, tail.port_a) annotation(
        Line(points = {{90, 50}, {100, 50}, {100, -20}}, color = {0, 127, 255}));
      connect(tail.port_a, valveTailOut.port_a) annotation(
        Line(points = {{100, -20}, {100, 20}, {110, 20}}, color = {0, 127, 255}));
      connect(pressurizedAir.ports[1], valveHeadIn.port_a) annotation(
        Line(points = {{-170, 80}, {-140, 80}, {-140, 50}, {-130, 50}}, color = {0, 127, 255}));
      connect(pressurizedAir.ports[2], valveBodyIn.port_a) annotation(
        Line(points = {{-170, 80}, {-40, 80}, {-40, 50}, {-30, 50}}, color = {0, 127, 255}));
      connect(pressurizedAir.ports[3], valveTailIn.port_a) annotation(
        Line(points = {{-170, 80}, {60, 80}, {60, 50}, {70, 50}}, color = {0, 127, 255}));
      connect(valveHeadOut.port_b, ambient.ports[1]) annotation(
        Line(points = {{-70, 18}, {-60, 18}, {-60, 100}, {150, 100}}, color = {0, 127, 255}));
      connect(valveBodyOut.port_b, ambient.ports[2]) annotation(
        Line(points = {{30, 20}, {40, 20}, {40, 100}, {150, 100}}, color = {0, 127, 255}));
      connect(valveTailOut.port_b, ambient.ports[3]) annotation(
        Line(points = {{130, 20}, {140, 20}, {140, 100}, {150, 100}}, color = {0, 127, 255}));
      connect(control.pin_n, ground.p) annotation(
        Line(points = {{-140, 130}, {-160, 130}}, color = {255, 0, 0}, thickness = 0.5));
      connect(control.pin_p[1], valveHeadIn.pin_p) annotation(
        Line(points = {{-120, 130}, {-110, 130}, {-110, 100}, {-136, 100}, {-136, 56}, {-130, 56}}, color = {255, 0, 0}));
      connect(control.pin_p[2], valveHeadOut.pin_p) annotation(
        Line(points = {{-120, 130}, {-96, 130}, {-96, 24}, {-90, 24}}, color = {255, 0, 0}));
      connect(control.pin_p[3], valveBodyIn.pin_p) annotation(
        Line(points = {{-120, 130}, {-36, 130}, {-36, 56}, {-30, 56}}, color = {255, 0, 0}));
      connect(control.pin_p[4], valveBodyOut.pin_p) annotation(
        Line(points = {{-120, 130}, {4, 130}, {4, 26}, {10, 26}}, color = {255, 0, 0}));
      connect(control.pin_p[5], valveTailIn.pin_p) annotation(
        Line(points = {{-120, 130}, {64, 130}, {64, 56}, {70, 56}}, color = {255, 0, 0}));
      connect(control.pin_p[6], valveTailOut.pin_p) annotation(
        Line(points = {{-120, 130}, {104, 130}, {104, 26}, {110, 26}}, color = {255, 0, 0}));
      connect(valveHeadIn.pin_n, ground[1].p) annotation(
        Line(points = {{-110, 56}, {-100, 56}, {-100, 110}, {-160, 110}, {-160, 130}}, color = {255, 0, 0}));
      connect(valveHeadOut.pin_n, ground[2].p) annotation(
        Line(points = {{-70, 24}, {-64, 24}, {-64, 110}, {-160, 110}, {-160, 130}}, color = {255, 0, 0}));
      connect(valveBodyIn.pin_n, ground[3].p) annotation(
        Line(points = {{-10, 56}, {0, 56}, {0, 110}, {-160, 110}, {-160, 130}}, color = {255, 0, 0}));
      connect(valveBodyOut.pin_n, ground[4].p) annotation(
        Line(points = {{30, 26}, {36, 26}, {36, 110}, {-160, 110}, {-160, 130}}, color = {255, 0, 0}));
      connect(valveTailIn.pin_n, ground[5].p) annotation(
        Line(points = {{90, 56}, {100, 56}, {100, 110}, {-160, 110}, {-160, 130}}, color = {255, 0, 0}));
      connect(valveTailOut.pin_n, ground[6].p) annotation(
        Line(points = {{130, 26}, {136, 26}, {136, 110}, {-160, 110}, {-160, 130}}, color = {255, 0, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-200, 140}, {180, -60}})),
        experiment(StartTime = 0, StopTime = 671, Tolerance = 1e-06, Interval = 0.1));
    end WormSystem2;

    model WormStatic
      replaceable package Medium = Modelica.Media.Air.DryAirNasa;
      //
      // Parameters
      //
      parameter Modelica.Units.SI.Mass m = 1;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 100000 "Nominal pressure difference across the valve";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.00237 "Nominal mass flow rate across the valve at dp_nominal";
      //
      // Components
      //
      Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
        Placement(visible = true, transformation(origin = {-196, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
        Placement(visible = true, transformation(origin = {-70, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 3, p = Medium.p_default*1.1) annotation(
        Placement(visible = true, transformation(origin = {-226, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.CylindricalActuators.LinearActuatorStaticLS bodyLinear(redeclare package Medium = Medium, l_0 = 0.127, r_0 = 0.0375) annotation(
        Placement(visible = true, transformation(origin = {1.77636e-15, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.SphericalActuators.RadialPressure tailRadial(redeclare package Medium = Medium, gamma_max = 2, m = 1, p_start = 2*system.p_start) annotation(
        Placement(visible = true, transformation(origin = {160, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3) annotation(
        Placement(visible = true, transformation(origin = {236, 62}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      WormBot.Components.SphericalActuators.RadialPressure headRadial(redeclare package Medium = Medium, gamma_max = 2, m = 1) annotation(
        Placement(visible = true, transformation(origin = {-140, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.ValveControl.LUTControl lUTControl annotation(
        Placement(visible = true, transformation(origin = {-210, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-270, 110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      WormBot.Components.Valves.SolenoidValve headValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-166, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve headValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-98, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {36, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {198, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
        Line(points = {{-20, -8}, {-120, -8}}, color = {0, 127, 0}));
      connect(headRadial.flange_a, mass.flange_a) annotation(
        Line(points = {{-160, -8}, {-186, -8}}, color = {0, 127, 0}));
      connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
        Line(points = {{20, -8}, {140, -8}}, color = {0, 127, 0}));
      connect(lUTControl.pin_n, ground.p) annotation(
        Line(points = {{-220, 120.2}, {-260, 120.2}, {-260, 110.2}}, color = {255, 0, 0}, thickness = 0.5));
      connect(headValveIn.port_b, headRadial.port_a) annotation(
        Line(points = {{-156, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(headValveIn.port_a, pressurizedAir.ports[1]) annotation(
        Line(points = {{-176, 34}, {-196, 34}, {-196, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(headValveIn.pin_n, ground[1].p) annotation(
        Line(points = {{-156, 44}, {-136, 44}, {-136, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[2], headValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {-116, 120}, {-116, 44}, {-108, 44}}, color = {255, 0, 0}));
      connect(headValveOut.pin_n, ground[2].p) annotation(
        Line(points = {{-88, 44}, {-76, 44}, {-76, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(pressurizedAir.ports[2], bodyValveIn.port_a) annotation(
        Line(points = {{-216, 72}, {-60, 72}, {-60, 34}, {-42, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_a, pressurizedAir.ports[3]) annotation(
        Line(points = {{120, 34}, {100, 34}, {100, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(bodyValveIn.port_b, bodyLinear.port_a) annotation(
        Line(points = {{-22, 34}, {0, 34}, {0, 12}}, color = {0, 127, 255}));
      connect(bodyLinear.port_a, bodyValveOut.port_a) annotation(
        Line(points = {{1.77636e-15, 12}, {1.77636e-15, 34}, {26, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_b, tailRadial.port_a) annotation(
        Line(points = {{140, 34}, {160, 34}, {160, 12}}, color = {0, 127, 255}));
      connect(tailRadial.port_a, tailValveOut.port_a) annotation(
        Line(points = {{160, 12}, {160, 34}, {188, 34}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[3], bodyValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-52, 120}, {-52, 44}, {-42, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[4], bodyValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {16, 120}, {16, 44}, {26, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[5], tailValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {108, 120}, {108, 44}, {120, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[6], tailValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {178, 120}, {178, 44}, {188, 44}}, color = {255, 0, 0}));
      connect(bodyValveIn.pin_n, ground[3].p) annotation(
        Line(points = {{-22, 44}, {-12, 44}, {-12, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(bodyValveOut.pin_n, ground[4].p) annotation(
        Line(points = {{46, 44}, {70, 44}, {70, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveIn.pin_n, ground[5].p) annotation(
        Line(points = {{140, 44}, {152, 44}, {152, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveOut.pin_n, ground[6].p) annotation(
        Line(points = {{208, 44}, {216, 44}, {216, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(headValveOut.port_b, ambient.ports[1]) annotation(
        Line(points = {{-88, 34}, {-68, 34}, {-68, 62}, {226, 62}}, color = {0, 127, 255}));
      connect(bodyValveOut.port_b, ambient.ports[2]) annotation(
        Line(points = {{46, 34}, {84, 34}, {84, 62}, {226, 62}}, color = {0, 127, 255}));
      connect(tailValveOut.port_b, ambient.ports[3]) annotation(
        Line(points = {{208, 34}, {222, 34}, {222, 62}, {226, 62}}, color = {0, 127, 255}));
      connect(headValveOut.port_a, headRadial.port_a) annotation(
        Line(points = {{-108, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[1], headValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-186, 120}, {-186, 44}, {-176, 44}}, color = {255, 0, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 75, Tolerance = 1e-06, Interval = 0.01),
        Diagram(coordinateSystem(extent = {{-280, 180}, {260, -40}})),
        Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
    end WormStatic;

    model WormDynamic
      replaceable package Medium = Modelica.Media.Air.DryAirNasa;
      //
      // Parameters
      //
      parameter Modelica.Units.SI.Mass m = 1;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 100000 "Nominal pressure difference across the valve";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.00237 "Nominal mass flow rate across the valve at dp_nominal";
      //
      // Components
      //
      Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
        Placement(visible = true, transformation(origin = {-196, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
        Placement(visible = true, transformation(origin = {-70, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 3, p = Medium.p_default*1.5) annotation(
        Placement(visible = true, transformation(origin = {-226, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.CylindricalActuators.LinearActuatorDynamic bodyLinear(redeclare package Medium = Medium, l_0 = 0.150, r_0 = 0.0375, t = 0.002) annotation(
        Placement(visible = true, transformation(origin = {1.77636e-15, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.SphericalActuators.RadialPressure tailRadial(redeclare package Medium = Medium, gamma_max = 2, p_start = 2*system.p_start) annotation(
        Placement(visible = true, transformation(origin = {160, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3) annotation(
        Placement(visible = true, transformation(origin = {236, 62}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      WormBot.Components.SphericalActuators.RadialPressure headRadial(redeclare package Medium = Medium, gamma_max = 2) annotation(
        Placement(visible = true, transformation(origin = {-140, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.ValveControl.LUTControl lUTControl annotation(
        Placement(visible = true, transformation(origin = {-210, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-270, 110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      WormBot.Components.Valves.SolenoidValve headValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-166, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve headValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-98, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {36, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {198, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
        Line(points = {{-20, -8}, {-120, -8}}, color = {0, 127, 0}));
      connect(headRadial.flange_a, mass.flange_a) annotation(
        Line(points = {{-160, -8}, {-186, -8}}, color = {0, 127, 0}));
      connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
        Line(points = {{20, -8}, {140, -8}}, color = {0, 127, 0}));
      connect(lUTControl.pin_n, ground.p) annotation(
        Line(points = {{-220, 120.2}, {-260, 120.2}, {-260, 110.2}}, color = {255, 0, 0}, thickness = 0.5));
      connect(headValveIn.port_b, headRadial.port_a) annotation(
        Line(points = {{-156, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(headValveIn.port_a, pressurizedAir.ports[1]) annotation(
        Line(points = {{-176, 34}, {-196, 34}, {-196, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(headValveIn.pin_n, ground[1].p) annotation(
        Line(points = {{-156, 44}, {-136, 44}, {-136, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[2], headValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {-116, 120}, {-116, 44}, {-108, 44}}, color = {255, 0, 0}));
      connect(headValveOut.pin_n, ground[2].p) annotation(
        Line(points = {{-88, 44}, {-76, 44}, {-76, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(pressurizedAir.ports[2], bodyValveIn.port_a) annotation(
        Line(points = {{-216, 72}, {-60, 72}, {-60, 34}, {-42, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_a, pressurizedAir.ports[3]) annotation(
        Line(points = {{120, 34}, {100, 34}, {100, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(bodyValveIn.port_b, bodyLinear.port_a) annotation(
        Line(points = {{-22, 34}, {0, 34}, {0, 12}}, color = {0, 127, 255}));
      connect(bodyLinear.port_a, bodyValveOut.port_a) annotation(
        Line(points = {{1.77636e-15, 12}, {1.77636e-15, 34}, {26, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_b, tailRadial.port_a) annotation(
        Line(points = {{140, 34}, {160, 34}, {160, 12}}, color = {0, 127, 255}));
      connect(tailRadial.port_a, tailValveOut.port_a) annotation(
        Line(points = {{160, 12}, {160, 34}, {188, 34}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[3], bodyValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-52, 120}, {-52, 44}, {-42, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[4], bodyValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {16, 120}, {16, 44}, {26, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[5], tailValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {108, 120}, {108, 44}, {120, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[6], tailValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {178, 120}, {178, 44}, {188, 44}}, color = {255, 0, 0}));
      connect(bodyValveIn.pin_n, ground[3].p) annotation(
        Line(points = {{-22, 44}, {-12, 44}, {-12, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(bodyValveOut.pin_n, ground[4].p) annotation(
        Line(points = {{46, 44}, {70, 44}, {70, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveIn.pin_n, ground[5].p) annotation(
        Line(points = {{140, 44}, {152, 44}, {152, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveOut.pin_n, ground[6].p) annotation(
        Line(points = {{208, 44}, {216, 44}, {216, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(headValveOut.port_b, ambient.ports[1]) annotation(
        Line(points = {{-88, 34}, {-68, 34}, {-68, 62}, {226, 62}}, color = {0, 127, 255}));
      connect(bodyValveOut.port_b, ambient.ports[2]) annotation(
        Line(points = {{46, 34}, {84, 34}, {84, 62}, {226, 62}}, color = {0, 127, 255}));
      connect(tailValveOut.port_b, ambient.ports[3]) annotation(
        Line(points = {{208, 34}, {222, 34}, {222, 62}, {226, 62}}, color = {0, 127, 255}));
      connect(headValveOut.port_a, headRadial.port_a) annotation(
        Line(points = {{-108, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[1], headValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-186, 120}, {-186, 44}, {-176, 44}}, color = {255, 0, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 350, Tolerance = 1e-06, Interval = 0.01),
        Diagram(coordinateSystem(extent = {{-280, 180}, {260, -40}})),
        Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
    end WormDynamic;

    model WormStaticLS
      replaceable package Medium = Modelica.Media.Air.DryAirNasa;
      //
      // Parameters
      //
      parameter Modelica.Units.SI.Mass m = 1;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 100000 "Nominal pressure difference across the valve";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.00237 "Nominal mass flow rate across the valve at dp_nominal";
      //
      // Components
      //
      Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
        Placement(visible = true, transformation(origin = {-196, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
        Placement(visible = true, transformation(origin = {-70, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 3, p = Medium.p_default*1.5) annotation(
        Placement(visible = true, transformation(origin = {-226, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.CylindricalActuators.LinearActuatorStatic bodyLinear(redeclare package Medium = Medium, L = 0.127, R = 0.0375, t = 0.002) annotation(
        Placement(visible = true, transformation(origin = {1.77636e-15, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.SphericalActuators.RadialPressure tailRadial(redeclare package Medium = Medium, gamma_max = 2, m = 1, p_start = 2*system.p_start) annotation(
        Placement(visible = true, transformation(origin = {160, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3) annotation(
        Placement(visible = true, transformation(origin = {236, 62}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      WormBot.Components.SphericalActuators.RadialPressure headRadial(redeclare package Medium = Medium, gamma_max = 2, m = 1) annotation(
        Placement(visible = true, transformation(origin = {-140, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.ValveControl.LUTControl lUTControl annotation(
        Placement(visible = true, transformation(origin = {-210, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-270, 110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      WormBot.Components.Valves.SolenoidValve headValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-166, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve headValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-98, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {36, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {198, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
        Line(points = {{-20, -8}, {-120, -8}}, color = {0, 127, 0}));
      connect(headRadial.flange_a, mass.flange_a) annotation(
        Line(points = {{-160, -8}, {-186, -8}}, color = {0, 127, 0}));
      connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
        Line(points = {{20, -8}, {140, -8}}, color = {0, 127, 0}));
      connect(lUTControl.pin_n, ground.p) annotation(
        Line(points = {{-220, 120.2}, {-260, 120.2}, {-260, 110.2}}, color = {255, 0, 0}, thickness = 0.5));
      connect(headValveIn.port_b, headRadial.port_a) annotation(
        Line(points = {{-156, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(headValveIn.port_a, pressurizedAir.ports[1]) annotation(
        Line(points = {{-176, 34}, {-196, 34}, {-196, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(headValveIn.pin_n, ground[1].p) annotation(
        Line(points = {{-156, 44}, {-136, 44}, {-136, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[2], headValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {-116, 120}, {-116, 44}, {-108, 44}}, color = {255, 0, 0}));
      connect(headValveOut.pin_n, ground[2].p) annotation(
        Line(points = {{-88, 44}, {-76, 44}, {-76, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(pressurizedAir.ports[2], bodyValveIn.port_a) annotation(
        Line(points = {{-216, 72}, {-60, 72}, {-60, 34}, {-42, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_a, pressurizedAir.ports[3]) annotation(
        Line(points = {{120, 34}, {100, 34}, {100, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(bodyValveIn.port_b, bodyLinear.port_a) annotation(
        Line(points = {{-22, 34}, {0, 34}, {0, 12}}, color = {0, 127, 255}));
      connect(bodyLinear.port_a, bodyValveOut.port_a) annotation(
        Line(points = {{1.77636e-15, 12}, {1.77636e-15, 34}, {26, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_b, tailRadial.port_a) annotation(
        Line(points = {{140, 34}, {160, 34}, {160, 12}}, color = {0, 127, 255}));
      connect(tailRadial.port_a, tailValveOut.port_a) annotation(
        Line(points = {{160, 12}, {160, 34}, {188, 34}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[3], bodyValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-52, 120}, {-52, 44}, {-42, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[4], bodyValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {16, 120}, {16, 44}, {26, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[5], tailValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {108, 120}, {108, 44}, {120, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[6], tailValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {178, 120}, {178, 44}, {188, 44}}, color = {255, 0, 0}));
      connect(bodyValveIn.pin_n, ground[3].p) annotation(
        Line(points = {{-22, 44}, {-12, 44}, {-12, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(bodyValveOut.pin_n, ground[4].p) annotation(
        Line(points = {{46, 44}, {70, 44}, {70, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveIn.pin_n, ground[5].p) annotation(
        Line(points = {{140, 44}, {152, 44}, {152, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveOut.pin_n, ground[6].p) annotation(
        Line(points = {{208, 44}, {216, 44}, {216, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(headValveOut.port_b, ambient.ports[1]) annotation(
        Line(points = {{-88, 34}, {-68, 34}, {-68, 62}, {226, 62}}, color = {0, 127, 255}));
      connect(bodyValveOut.port_b, ambient.ports[2]) annotation(
        Line(points = {{46, 34}, {84, 34}, {84, 62}, {226, 62}}, color = {0, 127, 255}));
      connect(tailValveOut.port_b, ambient.ports[3]) annotation(
        Line(points = {{208, 34}, {222, 34}, {222, 62}, {226, 62}}, color = {0, 127, 255}));
      connect(headValveOut.port_a, headRadial.port_a) annotation(
        Line(points = {{-108, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[1], headValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-186, 120}, {-186, 44}, {-176, 44}}, color = {255, 0, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 350, Tolerance = 1e-06, Interval = 0.01),
        Diagram(coordinateSystem(extent = {{-280, 180}, {260, -40}})),
        Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
    end WormStaticLS;

    model WormStaticMflow
      replaceable package Medium = Modelica.Media.Air.DryAirNasa;
      //
      // Parameters
      //
      parameter Modelica.Units.SI.Mass m = 1;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 100000 "Nominal pressure difference across the valve";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.00237 "Nominal mass flow rate across the valve at dp_nominal";
      //
      // Components
      //
      Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
        Placement(visible = true, transformation(origin = {-196, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
        Placement(visible = true, transformation(origin = {-70, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.MassFlowSource_T pressurizedAir(redeclare package Medium = Medium, m_flow = 0.017, nPorts = 4) annotation(
        Placement(visible = true, transformation(origin = {-226, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.CylindricalActuators.LinearActuatorStatic bodyLinear(redeclare package Medium = Medium, L = 0.127, R = 0.0375, t = 0.002) annotation(
        Placement(visible = true, transformation(origin = {1.77636e-15, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.SphericalActuators.RadialPressure tailRadial(redeclare package Medium = Medium, gamma_max = 2, m = 1, p_start = 2*system.p_start) annotation(
        Placement(visible = true, transformation(origin = {160, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 4) annotation(
        Placement(visible = true, transformation(origin = {352, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      WormBot.Components.SphericalActuators.RadialPressure headRadial(redeclare package Medium = Medium, gamma_max = 2, m = 1) annotation(
        Placement(visible = true, transformation(origin = {-140, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.ValveControl.LUTControl lUTControl annotation(
        Placement(visible = true, transformation(origin = {-210, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-270, 110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      WormBot.Components.Valves.SolenoidValve headValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-166, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve headValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-98, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {36, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {198, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Valves.ValveLinear valveLinear(redeclare package Medium = Medium, dp_nominal = 49999.99999999999, m_flow_nominal = 0.01) annotation(
        Placement(visible = true, transformation(origin = {316, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = 0.01, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nPorts = 3, use_portsData = false) annotation(
        Placement(visible = true, transformation(origin = {250, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {280, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Logical.OnOffController onOffController(bandwidth = 50000) annotation(
        Placement(visible = true, transformation(origin = {320, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression p_lim(y = 175000) annotation(
        Placement(visible = true, transformation(origin = {280, 180}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Logical.Switch switch1 annotation(
        Placement(visible = true, transformation(origin = {370, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression(y = 0) annotation(
        Placement(visible = true, transformation(origin = {330, 180}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression1(y = 1) annotation(
        Placement(visible = true, transformation(origin = {330, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
        Line(points = {{-20, -8}, {-120, -8}}, color = {0, 127, 0}));
      connect(headRadial.flange_a, mass.flange_a) annotation(
        Line(points = {{-160, -8}, {-186, -8}}, color = {0, 127, 0}));
      connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
        Line(points = {{20, -8}, {140, -8}}, color = {0, 127, 0}));
      connect(lUTControl.pin_n, ground.p) annotation(
        Line(points = {{-220, 120.2}, {-260, 120.2}, {-260, 110.2}}, color = {255, 0, 0}, thickness = 0.5));
      connect(headValveIn.port_b, headRadial.port_a) annotation(
        Line(points = {{-156, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(headValveIn.port_a, pressurizedAir.ports[1]) annotation(
        Line(points = {{-176, 34}, {-196, 34}, {-196, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(headValveIn.pin_n, ground[1].p) annotation(
        Line(points = {{-156, 44}, {-136, 44}, {-136, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[2], headValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {-116, 120}, {-116, 44}, {-108, 44}}, color = {255, 0, 0}));
      connect(headValveOut.pin_n, ground[2].p) annotation(
        Line(points = {{-88, 44}, {-76, 44}, {-76, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(pressurizedAir.ports[2], bodyValveIn.port_a) annotation(
        Line(points = {{-216, 72}, {-60, 72}, {-60, 34}, {-42, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_a, pressurizedAir.ports[3]) annotation(
        Line(points = {{120, 34}, {100, 34}, {100, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(bodyValveIn.port_b, bodyLinear.port_a) annotation(
        Line(points = {{-22, 34}, {0, 34}, {0, 12}}, color = {0, 127, 255}));
      connect(bodyLinear.port_a, bodyValveOut.port_a) annotation(
        Line(points = {{1.77636e-15, 12}, {1.77636e-15, 34}, {26, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_b, tailRadial.port_a) annotation(
        Line(points = {{140, 34}, {160, 34}, {160, 12}}, color = {0, 127, 255}));
      connect(tailRadial.port_a, tailValveOut.port_a) annotation(
        Line(points = {{160, 12}, {160, 34}, {188, 34}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[3], bodyValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-52, 120}, {-52, 44}, {-42, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[4], bodyValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {16, 120}, {16, 44}, {26, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[5], tailValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {108, 120}, {108, 44}, {120, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[6], tailValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {178, 120}, {178, 44}, {188, 44}}, color = {255, 0, 0}));
      connect(bodyValveIn.pin_n, ground[3].p) annotation(
        Line(points = {{-22, 44}, {-12, 44}, {-12, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(bodyValveOut.pin_n, ground[4].p) annotation(
        Line(points = {{46, 44}, {70, 44}, {70, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveIn.pin_n, ground[5].p) annotation(
        Line(points = {{140, 44}, {152, 44}, {152, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveOut.pin_n, ground[6].p) annotation(
        Line(points = {{208, 44}, {216, 44}, {216, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(headValveOut.port_b, ambient.ports[1]) annotation(
        Line(points = {{-88, 34}, {-68, 34}, {-68, 60}, {342, 60}}, color = {0, 127, 255}));
      connect(bodyValveOut.port_b, ambient.ports[2]) annotation(
        Line(points = {{46, 34}, {84, 34}, {84, 60}, {342, 60}}, color = {0, 127, 255}));
      connect(tailValveOut.port_b, ambient.ports[3]) annotation(
        Line(points = {{208, 34}, {222, 34}, {222, 60}, {342, 60}}, color = {0, 127, 255}));
      connect(headValveOut.port_a, headRadial.port_a) annotation(
        Line(points = {{-108, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[1], headValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-186, 120}, {-186, 44}, {-176, 44}}, color = {255, 0, 0}));
      connect(volume.ports[1], valveLinear.port_a) annotation(
        Line(points = {{250, 100}, {250, 90}, {306, 90}}, color = {0, 127, 255}));
      connect(pressurizedAir.ports[4], volume.ports[2]) annotation(
        Line(points = {{-216, 72}, {250, 72}, {250, 100}}, color = {0, 127, 255}));
      connect(valveLinear.port_b, ambient.ports[4]) annotation(
        Line(points = {{326, 90}, {334, 90}, {334, 60}, {342, 60}}, color = {0, 127, 255}));
      connect(volume.ports[3], pressure.port) annotation(
        Line(points = {{250, 100}, {250, 90}, {280, 90}, {280, 140}}, color = {0, 127, 255}));
      connect(pressure.p, onOffController.u) annotation(
        Line(points = {{292, 150}, {300, 150}, {300, 144}, {308, 144}}, color = {0, 0, 127}));
      connect(p_lim.y, onOffController.reference) annotation(
        Line(points = {{292, 180}, {300, 180}, {300, 156}, {308, 156}}, color = {0, 0, 127}));
      connect(onOffController.y, switch1.u2) annotation(
        Line(points = {{332, 150}, {358, 150}}, color = {255, 0, 255}));
      connect(switch1.y, valveLinear.opening) annotation(
        Line(points = {{382, 150}, {400, 150}, {400, 110}, {316, 110}, {316, 98}}, color = {0, 0, 127}));
      connect(realExpression1.y, switch1.u3) annotation(
        Line(points = {{342, 120}, {346, 120}, {346, 142}, {358, 142}}, color = {0, 0, 127}));
      connect(realExpression.y, switch1.u1) annotation(
        Line(points = {{342, 180}, {348, 180}, {348, 158}, {358, 158}}, color = {0, 0, 127}));
      annotation(
        experiment(StartTime = 0, StopTime = 75, Tolerance = 1e-06, Interval = 0.01),
        Diagram(coordinateSystem(extent = {{-280, 200}, {400, -40}})),
        Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
    end WormStaticMflow;

    model WormStaticMflow1
      replaceable package Medium = Modelica.Media.Air.DryAirNasa;
      //
      // Parameters
      //
      parameter Modelica.Units.SI.Mass m = 1;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 100000 "Nominal pressure difference across the valve";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.0015 "Nominal mass flow rate across the valve at dp_nominal";
      parameter Modelica.Units.SI.MassFlowRate m_flow_regulator = 0.000015 "Mass flow through the regulator";
      //
      // Components
      //
      Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
        Placement(visible = true, transformation(origin = {-196, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, p_ambient = 100567) annotation(
        Placement(visible = true, transformation(origin = {-70, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.MassFlowSource_T pressurizedAir[3](redeclare package Medium = Medium, each nPorts = 1, each use_m_flow_in = true) annotation(
        Placement(visible = true, transformation(origin = {-226, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.CylindricalActuators.LinearActuatorStaticLS bodyLinear(redeclare package Medium = Medium, l_0 = 0.146, p_ref = system.p_ambient, r_0 = 0.0375) annotation(
        Placement(visible = true, transformation(origin = {1.77636e-15, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.RadialActuators.RadialPressureLS tailRadial(redeclare package Medium = Medium, m = 1, p_start = system.p_ambient) annotation(
        Placement(visible = true, transformation(origin = {160, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3, p = system.p_ambient) annotation(
        Placement(visible = true, transformation(origin = {250, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      WormBot.Components.RadialActuators.RadialPressureLS headRadial(redeclare package Medium = Medium, m = 1, p_ambient = system.p_ambient) annotation(
        Placement(visible = true, transformation(origin = {-140, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.ValveControl.LUTControl lUTControl annotation(
        Placement(visible = true, transformation(origin = {-210, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-270, 110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      WormBot.Components.Valves.SolenoidValve headValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-166, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve headValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-98, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {36, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {198, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Gain gain[3](each k = m_flow_regulator) annotation(
        Placement(visible = true, transformation(origin = {-280, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
        Line(points = {{-20, -8}, {-120, -8}}, color = {0, 127, 0}));
      connect(headRadial.flange_a, mass.flange_a) annotation(
        Line(points = {{-160, -8}, {-186, -8}}, color = {0, 127, 0}));
      connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
        Line(points = {{20, -8}, {140, -8}}, color = {0, 127, 0}));
      connect(lUTControl.pin_n, ground.p) annotation(
        Line(points = {{-220, 120.2}, {-260, 120.2}, {-260, 110.2}}, color = {255, 0, 0}, thickness = 0.5));
      connect(headValveIn.port_b, headRadial.port_a) annotation(
        Line(points = {{-156, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(headValveIn.port_a, pressurizedAir[1].ports[1]) annotation(
        Line(points = {{-176, 34}, {-196, 34}, {-196, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(headValveIn.pin_n, ground[1].p) annotation(
        Line(points = {{-156, 44}, {-136, 44}, {-136, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[2], headValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {-116, 120}, {-116, 44}, {-108, 44}}, color = {255, 0, 0}));
      connect(headValveOut.pin_n, ground[2].p) annotation(
        Line(points = {{-88, 44}, {-76, 44}, {-76, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(pressurizedAir[2].ports[1], bodyValveIn.port_a) annotation(
        Line(points = {{-216, 72}, {-60, 72}, {-60, 34}, {-42, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_a, pressurizedAir[3].ports[1]) annotation(
        Line(points = {{120, 34}, {100, 34}, {100, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(bodyValveIn.port_b, bodyLinear.port_a) annotation(
        Line(points = {{-22, 34}, {0, 34}, {0, 12}}, color = {0, 127, 255}));
      connect(bodyLinear.port_a, bodyValveOut.port_a) annotation(
        Line(points = {{1.77636e-15, 12}, {1.77636e-15, 34}, {26, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_b, tailRadial.port_a) annotation(
        Line(points = {{140, 34}, {160, 34}, {160, 12}}, color = {0, 127, 255}));
      connect(tailRadial.port_a, tailValveOut.port_a) annotation(
        Line(points = {{160, 12}, {160, 34}, {188, 34}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[3], bodyValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-52, 120}, {-52, 44}, {-42, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[4], bodyValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {16, 120}, {16, 44}, {26, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[5], tailValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {108, 120}, {108, 44}, {120, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[6], tailValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {178, 120}, {178, 44}, {188, 44}}, color = {255, 0, 0}));
      connect(bodyValveIn.pin_n, ground[3].p) annotation(
        Line(points = {{-22, 44}, {-12, 44}, {-12, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(bodyValveOut.pin_n, ground[4].p) annotation(
        Line(points = {{46, 44}, {70, 44}, {70, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveIn.pin_n, ground[5].p) annotation(
        Line(points = {{140, 44}, {152, 44}, {152, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveOut.pin_n, ground[6].p) annotation(
        Line(points = {{208, 44}, {216, 44}, {216, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(headValveOut.port_b, ambient.ports[1]) annotation(
        Line(points = {{-88, 34}, {-68, 34}, {-68, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(bodyValveOut.port_b, ambient.ports[2]) annotation(
        Line(points = {{46, 34}, {84, 34}, {84, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(tailValveOut.port_b, ambient.ports[3]) annotation(
        Line(points = {{208, 34}, {222, 34}, {222, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(headValveOut.port_a, headRadial.port_a) annotation(
        Line(points = {{-108, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[1], headValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-186, 120}, {-186, 44}, {-176, 44}}, color = {255, 0, 0}));
      connect(gain.y, pressurizedAir.m_flow_in) annotation(
        Line(points = {{-268, 72}, {-260, 72}, {-260, 80}, {-236, 80}}, color = {0, 0, 127}, thickness = 0.5));
      gain[1].u = lUTControl.controllerLUT.y[1];
      gain[2].u = lUTControl.controllerLUT.y[3];
      gain[3].u = lUTControl.controllerLUT.y[5];
      annotation(
        experiment(StartTime = 0, StopTime = 75, Tolerance = 1e-06, Interval = 0.01),
        Diagram(coordinateSystem(extent = {{-280, 200}, {400, -40}})),
        Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
    end WormStaticMflow1;

    model WormPressureIn
      replaceable package Medium = Modelica.Media.Air.DryAirNasa;
      //
      // Parameters
      //
      parameter Modelica.Units.SI.Mass m = 1;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 100000 "Nominal pressure difference across the valve";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.0015 "Nominal mass flow rate across the valve at dp_nominal";
      parameter Modelica.Units.SI.MassFlowRate m_flow_regulator = 0.000010 "Mass flow through the regulator";
      //
      // Components
      //
      Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
        Placement(visible = true, transformation(origin = {-196, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
        Placement(visible = true, transformation(origin = {-70, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.MassFlowSource_T pressurizedAir[3](redeclare package Medium = Medium, each nPorts = 1, each use_m_flow_in = true) annotation(
        Placement(visible = true, transformation(origin = {-226, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.CylindricalActuators.LinearActuatorStaticLS bodyLinear(redeclare package Medium = Medium, l_0 = 0.125, r_0 = 0.0375) annotation(
        Placement(visible = true, transformation(origin = {1.77636e-15, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.RadialActuators.RadialPressureContact tailRadial(redeclare package Medium = Medium, gamma_max = 2, m = 1, p_start = 2*system.p_start) annotation(
        Placement(visible = true, transformation(origin = {160, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3) annotation(
        Placement(visible = true, transformation(origin = {250, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      WormBot.Components.RadialActuators.RadialPressureContact headRadial(redeclare package Medium = Medium, gamma_max = 2, m = 1) annotation(
        Placement(visible = true, transformation(origin = {-140, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.ValveControl.LUTControl lUTControl annotation(
        Placement(visible = true, transformation(origin = {-210, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-270, 110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      WormBot.Components.Valves.SolenoidValve headValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-166, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve headValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-98, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {36, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {198, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Gain gain[3](each k = m_flow_regulator) annotation(
        Placement(visible = true, transformation(origin = {-280, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
        Line(points = {{-20, -8}, {-120, -8}}, color = {0, 127, 0}));
      connect(headRadial.flange_a, mass.flange_a) annotation(
        Line(points = {{-160, -8}, {-186, -8}}, color = {0, 127, 0}));
      connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
        Line(points = {{20, -8}, {140, -8}}, color = {0, 127, 0}));
      connect(lUTControl.pin_n, ground.p) annotation(
        Line(points = {{-220, 120.2}, {-260, 120.2}, {-260, 110.2}}, color = {255, 0, 0}, thickness = 0.5));
      connect(headValveIn.port_b, headRadial.port_a) annotation(
        Line(points = {{-156, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(headValveIn.port_a, pressurizedAir[1].ports[1]) annotation(
        Line(points = {{-176, 34}, {-196, 34}, {-196, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(headValveIn.pin_n, ground[1].p) annotation(
        Line(points = {{-156, 44}, {-136, 44}, {-136, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[2], headValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {-116, 120}, {-116, 44}, {-108, 44}}, color = {255, 0, 0}));
      connect(headValveOut.pin_n, ground[2].p) annotation(
        Line(points = {{-88, 44}, {-76, 44}, {-76, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(pressurizedAir[2].ports[1], bodyValveIn.port_a) annotation(
        Line(points = {{-216, 72}, {-60, 72}, {-60, 34}, {-42, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_a, pressurizedAir[3].ports[1]) annotation(
        Line(points = {{120, 34}, {100, 34}, {100, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(bodyValveIn.port_b, bodyLinear.port_a) annotation(
        Line(points = {{-22, 34}, {0, 34}, {0, 12}}, color = {0, 127, 255}));
      connect(bodyLinear.port_a, bodyValveOut.port_a) annotation(
        Line(points = {{1.77636e-15, 12}, {1.77636e-15, 34}, {26, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_b, tailRadial.port_a) annotation(
        Line(points = {{140, 34}, {160, 34}, {160, 12}}, color = {0, 127, 255}));
      connect(tailRadial.port_a, tailValveOut.port_a) annotation(
        Line(points = {{160, 12}, {160, 34}, {188, 34}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[3], bodyValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-52, 120}, {-52, 44}, {-42, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[4], bodyValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {16, 120}, {16, 44}, {26, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[5], tailValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {108, 120}, {108, 44}, {120, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[6], tailValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {178, 120}, {178, 44}, {188, 44}}, color = {255, 0, 0}));
      connect(bodyValveIn.pin_n, ground[3].p) annotation(
        Line(points = {{-22, 44}, {-12, 44}, {-12, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(bodyValveOut.pin_n, ground[4].p) annotation(
        Line(points = {{46, 44}, {70, 44}, {70, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveIn.pin_n, ground[5].p) annotation(
        Line(points = {{140, 44}, {152, 44}, {152, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveOut.pin_n, ground[6].p) annotation(
        Line(points = {{208, 44}, {216, 44}, {216, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(headValveOut.port_b, ambient.ports[1]) annotation(
        Line(points = {{-88, 34}, {-68, 34}, {-68, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(bodyValveOut.port_b, ambient.ports[2]) annotation(
        Line(points = {{46, 34}, {84, 34}, {84, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(tailValveOut.port_b, ambient.ports[3]) annotation(
        Line(points = {{208, 34}, {222, 34}, {222, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(headValveOut.port_a, headRadial.port_a) annotation(
        Line(points = {{-108, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[1], headValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-186, 120}, {-186, 44}, {-176, 44}}, color = {255, 0, 0}));
      connect(gain.y, pressurizedAir.m_flow_in) annotation(
        Line(points = {{-268, 72}, {-260, 72}, {-260, 80}, {-236, 80}}, color = {0, 0, 127}, thickness = 0.5));
      gain[1].u = lUTControl.controllerLUT.y[1];
      gain[2].u = lUTControl.controllerLUT.y[3];
      gain[3].u = lUTControl.controllerLUT.y[5];
      annotation(
        experiment(StartTime = 0, StopTime = 75, Tolerance = 1e-06, Interval = 0.01),
        Diagram(coordinateSystem(extent = {{-280, 200}, {400, -40}})),
        Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
    end WormPressureIn;

    model WormDynamicMflow1
      replaceable package Medium = Modelica.Media.Air.DryAirNasa;
      //
      // Parameters
      //
      parameter Modelica.Units.SI.Mass m = 1;
      parameter Modelica.Units.SI.Mass m_r = 1;
      parameter Modelica.Units.SI.Mass m_l = 1;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 100000 "Nominal pressure difference across the valve";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.0015 "Nominal mass flow rate across the valve at dp_nominal";
      parameter Modelica.Units.SI.MassFlowRate m_flow_regulator = 0.000015 "Mass flow through the regulator";
      //
      // Components
      //
      Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
        Placement(visible = true, transformation(origin = {-196, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, p_ambient = 100567) annotation(
        Placement(visible = true, transformation(origin = {-70, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.MassFlowSource_T pressurizedAir[3](redeclare package Medium = Medium, each nPorts = 1, each use_m_flow_in = true) annotation(
        Placement(visible = true, transformation(origin = {-226, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.CylindricalActuators.LinearActuatorDynamicDirection bodyLinear(redeclare package Medium = Medium, l_0 = 0.150, m_l = m_l, m_r = 0, p_ambient = system.p_ambient, r_0 = 0.0375) annotation(
        Placement(visible = true, transformation(origin = {1.77636e-15, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.RadialActuators.RadialPressureDynamic tailRadial(redeclare package Medium = Medium, m = m_r, mu_k_joint = 0.1, mu_s_joint = 0.15, p_ambient = system.p_ambient, p_start = 2*system.p_start) annotation(
        Placement(visible = true, transformation(origin = {160, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3, p = system.p_ambient) annotation(
        Placement(visible = true, transformation(origin = {250, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      WormBot.Components.RadialActuators.RadialPressureDynamic headRadial(redeclare package Medium = Medium, m = m_r, mu_k_joint = 0.1, mu_s_joint = 0.15, p_ambient = system.p_ambient) annotation(
        Placement(visible = true, transformation(origin = {-140, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.ValveControl.LUTControl lUTControl annotation(
        Placement(visible = true, transformation(origin = {-210, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-270, 110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      WormBot.Components.Valves.SolenoidValve headValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-166, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve headValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-98, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {36, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {198, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Gain gain[3](each k = m_flow_regulator) annotation(
        Placement(visible = true, transformation(origin = {-282, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
//
// Equations
//
//pressurizedAir[1].m_flow_in = lUTControl.controllerLUT.y[1] / m_flow_regulator;
//pressurizedAir[2].m_flow_in = lUTControl.controllerLUT.y[3] / m_flow_regulator;
//pressurizedAir[3].m_flow_in = lUTControl.controllerLUT.y[5] / m_flow_regulator;
//
// Connections
//
      connect(lUTControl.pin_n, ground.p) annotation(
        Line(points = {{-220, 120.2}, {-260, 120.2}, {-260, 110.2}}, color = {255, 0, 0}, thickness = 0.5));
      connect(headValveIn.port_b, headRadial.port_a) annotation(
        Line(points = {{-156, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(headValveIn.port_a, pressurizedAir[1].ports[1]) annotation(
        Line(points = {{-176, 34}, {-196, 34}, {-196, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(headValveIn.pin_n, ground[1].p) annotation(
        Line(points = {{-156, 44}, {-136, 44}, {-136, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[2], headValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {-116, 120}, {-116, 44}, {-108, 44}}, color = {255, 0, 0}));
      connect(headValveOut.pin_n, ground[2].p) annotation(
        Line(points = {{-88, 44}, {-76, 44}, {-76, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(pressurizedAir[2].ports[1], bodyValveIn.port_a) annotation(
        Line(points = {{-216, 72}, {-60, 72}, {-60, 34}, {-42, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_a, pressurizedAir[3].ports[1]) annotation(
        Line(points = {{120, 34}, {100, 34}, {100, 72}, {-216, 72}}, color = {0, 127, 255}));
      connect(bodyValveIn.port_b, bodyLinear.port_a) annotation(
        Line(points = {{-22, 34}, {0, 34}, {0, 12}}, color = {0, 127, 255}));
      connect(bodyLinear.port_a, bodyValveOut.port_a) annotation(
        Line(points = {{0, 12}, {0, 34}, {26, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_b, tailRadial.port_a) annotation(
        Line(points = {{140, 34}, {160, 34}, {160, 12}}, color = {0, 127, 255}));
      connect(tailRadial.port_a, tailValveOut.port_a) annotation(
        Line(points = {{160, 12}, {160, 34}, {188, 34}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[3], bodyValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-52, 120}, {-52, 44}, {-42, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[4], bodyValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {16, 120}, {16, 44}, {26, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[5], tailValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {108, 120}, {108, 44}, {120, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[6], tailValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {178, 120}, {178, 44}, {188, 44}}, color = {255, 0, 0}));
      connect(bodyValveIn.pin_n, ground[3].p) annotation(
        Line(points = {{-22, 44}, {-12, 44}, {-12, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(bodyValveOut.pin_n, ground[4].p) annotation(
        Line(points = {{46, 44}, {70, 44}, {70, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveIn.pin_n, ground[5].p) annotation(
        Line(points = {{140, 44}, {152, 44}, {152, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveOut.pin_n, ground[6].p) annotation(
        Line(points = {{208, 44}, {216, 44}, {216, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(headValveOut.port_b, ambient.ports[1]) annotation(
        Line(points = {{-88, 34}, {-68, 34}, {-68, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(bodyValveOut.port_b, ambient.ports[2]) annotation(
        Line(points = {{46, 34}, {84, 34}, {84, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(tailValveOut.port_b, ambient.ports[3]) annotation(
        Line(points = {{208, 34}, {222, 34}, {222, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(headValveOut.port_a, headRadial.port_a) annotation(
        Line(points = {{-108, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[1], headValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-186, 120}, {-186, 44}, {-176, 44}}, color = {255, 0, 0}));
      connect(gain.y, pressurizedAir.m_flow_in) annotation(
        Line(points = {{-271, 80}, {-236, 80}}, color = {0, 0, 127}, thickness = 0.5));
      gain[1].u = lUTControl.controllerLUT.y[1];
      gain[2].u = lUTControl.controllerLUT.y[3];
      gain[3].u = lUTControl.controllerLUT.y[5];
      connect(tailRadial.flange_a, bodyLinear.flange_b) annotation(
        Line(points = {{140, -8}, {20, -8}}, color = {0, 127, 0}));
      connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
        Line(points = {{-20, -8}, {-120, -8}}, color = {0, 127, 0}));
      connect(headRadial.flange_a, mass.flange_a) annotation(
        Line(points = {{-160, -8}, {-186, -8}}, color = {0, 127, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 75, Tolerance = 1e-06, Interval = 0.001),
        Diagram(coordinateSystem(extent = {{-300, 160}, {260, -40}})),
        Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
    end WormDynamicMflow1;

    model WormStaticComp
      replaceable package Medium = Modelica.Media.Air.DryAirNasa;
      //
      // Parameters
      //
      parameter Modelica.Units.SI.Mass m_head = 0.4782;
      parameter Modelica.Units.SI.Mass m_tail = 0.4942;
      parameter Modelica.Units.SI.Mass m_l = 0.1681;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 100000 "Nominal pressure difference across the valve";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.0015 "Nominal mass flow rate across the valve at dp_nominal"; // 0.0015
      parameter Modelica.Units.SI.MassFlowRate m_flow_regulator = 0.0000175 "Mass flow through the regulator";       // 0.000015 gives 32cm in 5 cycles
      //
      // Components
      //
      Modelica.Mechanics.Translational.Components.Mass mass(m = m_head) annotation(
        Placement(visible = true, transformation(origin = {-196, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, p_ambient = 100567) annotation(
        Placement(visible = true, transformation(origin = {-70, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.CylindricalActuators.LinearActuatorStaticLS bodyLinear(redeclare package Medium = Medium, l_0 = 0.146, p_ref = system.p_ambient, r_0 = 0.0375) annotation(
        Placement(visible = true, transformation(origin = {1.77636e-15, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.RadialActuators.RadialPressureLS tailRadial(redeclare package Medium = Medium, m = m_tail, p_start = system.p_ambient, r_joint = 0.095) annotation(
        Placement(visible = true, transformation(origin = {160, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3, p = system.p_ambient) annotation(
        Placement(visible = true, transformation(origin = {250, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      WormBot.Components.RadialActuators.RadialPressureLS headRadial(redeclare package Medium = Medium, m = m_head, p_ambient = system.p_ambient, r_joint = 0.095) annotation(
        Placement(visible = true, transformation(origin = {-140, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.ValveControl.LUTControl lUTControl annotation(
        Placement(visible = true, transformation(origin = {-210, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-270, 110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      WormBot.Components.Valves.SolenoidValve headValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-166, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve headValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-98, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {36, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {198, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.AirIn airIn(m_flow_regulator = m_flow_regulator) annotation(
        Placement(visible = true, transformation(origin = {-250, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Components.Mass mass1(m = m_tail) annotation(
        Placement(visible = true, transformation(origin = {212, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
        Line(points = {{-20, -8}, {-120, -8}}, color = {0, 127, 0}));
      connect(headRadial.flange_a, mass.flange_a) annotation(
        Line(points = {{-160, -8}, {-186, -8}}, color = {0, 127, 0}));
      connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
        Line(points = {{20, -8}, {140, -8}}, color = {0, 127, 0}));
      connect(lUTControl.pin_n, ground.p) annotation(
        Line(points = {{-220, 120.2}, {-260, 120.2}, {-260, 110.2}}, color = {255, 0, 0}, thickness = 0.5));
      connect(headValveIn.port_b, headRadial.port_a) annotation(
        Line(points = {{-156, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(headValveIn.pin_n, ground[1].p) annotation(
        Line(points = {{-156, 44}, {-136, 44}, {-136, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[2], headValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {-116, 120}, {-116, 44}, {-108, 44}}, color = {255, 0, 0}));
      connect(headValveOut.pin_n, ground[2].p) annotation(
        Line(points = {{-88, 44}, {-76, 44}, {-76, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(bodyValveIn.port_b, bodyLinear.port_a) annotation(
        Line(points = {{-22, 34}, {0, 34}, {0, 12}}, color = {0, 127, 255}));
      connect(bodyLinear.port_a, bodyValveOut.port_a) annotation(
        Line(points = {{1.77636e-15, 12}, {1.77636e-15, 34}, {26, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_b, tailRadial.port_a) annotation(
        Line(points = {{140, 34}, {160, 34}, {160, 12}}, color = {0, 127, 255}));
      connect(tailRadial.port_a, tailValveOut.port_a) annotation(
        Line(points = {{160, 12}, {160, 34}, {188, 34}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[3], bodyValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-52, 120}, {-52, 44}, {-42, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[4], bodyValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {16, 120}, {16, 44}, {26, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[5], tailValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {108, 120}, {108, 44}, {120, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[6], tailValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {178, 120}, {178, 44}, {188, 44}}, color = {255, 0, 0}));
      connect(bodyValveIn.pin_n, ground[3].p) annotation(
        Line(points = {{-22, 44}, {-12, 44}, {-12, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(bodyValveOut.pin_n, ground[4].p) annotation(
        Line(points = {{46, 44}, {70, 44}, {70, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveIn.pin_n, ground[5].p) annotation(
        Line(points = {{140, 44}, {152, 44}, {152, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveOut.pin_n, ground[6].p) annotation(
        Line(points = {{208, 44}, {216, 44}, {216, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(headValveOut.port_b, ambient.ports[1]) annotation(
        Line(points = {{-88, 34}, {-68, 34}, {-68, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(bodyValveOut.port_b, ambient.ports[2]) annotation(
        Line(points = {{46, 34}, {84, 34}, {84, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(tailValveOut.port_b, ambient.ports[3]) annotation(
        Line(points = {{208, 34}, {222, 34}, {222, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(headValveOut.port_a, headRadial.port_a) annotation(
        Line(points = {{-108, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[1], headValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-186, 120}, {-186, 44}, {-176, 44}}, color = {255, 0, 0}));
      airIn.onOff[1] = lUTControl.controllerLUT.y[1];
      airIn.onOff[2] = lUTControl.controllerLUT.y[3];
      airIn.onOff[3] = lUTControl.controllerLUT.y[5];
      connect(airIn.port_b[1], headValveIn.port_a) annotation(
        Line(points = {{-240, 70}, {-200, 70}, {-200, 34}, {-176, 34}}, color = {0, 127, 255}));
      connect(airIn.port_b[2], bodyValveIn.port_a) annotation(
        Line(points = {{-240, 70}, {-60, 70}, {-60, 34}, {-42, 34}}, color = {0, 127, 255}));
      connect(airIn.port_b[3], tailValveIn.port_a) annotation(
        Line(points = {{-240, 70}, {100, 70}, {100, 34}, {120, 34}}, color = {0, 127, 255}));
  connect(mass1.flange_b, tailRadial.flange_b) annotation(
        Line(points = {{202, -8}, {180, -8}}, color = {0, 127, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 75, Tolerance = 1e-06, Interval = 0.01),
        Diagram(coordinateSystem(extent = {{-280, 200}, {400, -40}})),
        Documentation(info = "<html><head></head><body>This is the final complete worm model. This model reads in the .txt file for the worm control, so the filepath to that file may need to be adjusted on difference devices. The linear actuator uses the least-squares polynomial fit model for the elongation</body></html>"));
    end WormStaticComp;

    model WormDynamicComp
      replaceable package Medium = Modelica.Media.Air.DryAirNasa;
      //
      // Parameters
      //
      parameter Modelica.Units.SI.Mass m_head = 0.4782;
      parameter Modelica.Units.SI.Mass m_tail = 0.4942;
      parameter Modelica.Units.SI.Mass m_l = 0.1681;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 100000 "Nominal pressure difference across the valve";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.0015 "Nominal mass flow rate across the valve at dp_nominal";
      parameter Modelica.Units.SI.MassFlowRate m_flow_regulator = 0.0000175 "Mass flow through the regulator";
      //
      // Components
      //
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, p_ambient = 100567) annotation(
        Placement(visible = true, transformation(origin = {-67, 143}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
      WormBot.Components.CylindricalActuators.LinearActuatorDynamicStrain bodyLinear(redeclare package Medium = Medium, b = 50, l_0 = 0.150, m_l = m_l, m_r = (m_head+m_tail)/2, p_ambient = system.p_ambient, r_0 = 0.0375, t = 0.005) annotation(
        Placement(visible = true, transformation(origin = {1.77636e-15, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.RadialActuators.RadialPressureDynamic tailRadial(redeclare package Medium = Medium, m = m_tail, mu_k_actuator = 0.6, mu_k_joint = 0.2, mu_s_actuator = 0.65, mu_s_joint = 0.25, p_ambient = system.p_ambient, p_start = 2*system.p_start) annotation(
        Placement(visible = true, transformation(origin = {160, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3, p = system.p_ambient) annotation(
        Placement(visible = true, transformation(origin = {252, 63}, extent = {{20, -19}, {-20, 19}}, rotation = 0)));
      WormBot.Components.RadialActuators.RadialPressureDynamic headRadial(redeclare package Medium = Medium, m = m_head, mu_k_actuator = 0.6, mu_k_joint = 0.2, mu_s_actuator = 0.65, mu_s_joint = 0.25, p_ambient = system.p_ambient) annotation(
        Placement(visible = true, transformation(origin = {-140, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.ValveControl.LUTControl lUTControl annotation(
        Placement(visible = true, transformation(origin = {-214, 130}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-275, 115}, extent = {{-15, -15}, {15, 15}}, rotation = -90)));
      WormBot.Components.Valves.SolenoidValve headValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-162, 42}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve headValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-96, 42}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {42, 42}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-32, 42}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {198, 42}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {130, 42}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      WormBot.Components.AirIn airIn(m_flow_regulator = m_flow_regulator) annotation(
        Placement(visible = true, transformation(origin = {-240, 72}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
  Modelica.Mechanics.Translational.Components.Mass headMass(m = m_head) annotation(
        Placement(visible = true, transformation(origin = {-196, -8}, extent = {{20, -20}, {-20, 20}}, rotation = 0)));
  Modelica.Mechanics.Translational.Components.Mass tailMass(m = m_tail) annotation(
        Placement(visible = true, transformation(origin = {220, -8}, extent = {{20, -20}, {-20, 20}}, rotation = 0)));
    equation
//
// Equations
//
//pressurizedAir[1].m_flow_in = lUTControl.controllerLUT.y[1] / m_flow_regulator;
//pressurizedAir[2].m_flow_in = lUTControl.controllerLUT.y[3] / m_flow_regulator;
//pressurizedAir[3].m_flow_in = lUTControl.controllerLUT.y[5] / m_flow_regulator;
//
// Connections
//
      connect(lUTControl.pin_n, ground.p) annotation(
        Line(points = {{-230, 130.32}, {-260, 130.32}, {-260, 115.32}}, color = {255, 0, 0}, thickness = 0.5));
      connect(headValveIn.port_b, headRadial.port_a) annotation(
        Line(points = {{-148, 36}, {-140, 36}, {-140, 12}}, color = {0, 127, 255}));
      connect(headValveIn.pin_n, ground[1].p) annotation(
        Line(points = {{-148, 50}, {-136, 50}, {-136, 100}, {-260, 100}, {-260, 115}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[2], headValveOut.pin_p) annotation(
        Line(points = {{-198, 130}, {-116, 130}, {-116, 50}, {-110, 50}}, color = {255, 0, 0}));
      connect(headValveOut.pin_n, ground[2].p) annotation(
        Line(points = {{-82, 50}, {-76, 50}, {-76, 100}, {-260, 100}, {-260, 115}}, color = {255, 0, 0}));
      connect(bodyValveIn.port_b, bodyLinear.port_a) annotation(
        Line(points = {{-18, 36}, {0, 36}, {0, 12}}, color = {0, 127, 255}));
      connect(bodyLinear.port_a, bodyValveOut.port_a) annotation(
        Line(points = {{0, 12}, {0, 36}, {28, 36}}, color = {0, 127, 255}));
      connect(tailValveIn.port_b, tailRadial.port_a) annotation(
        Line(points = {{144, 36}, {160, 36}, {160, 12}}, color = {0, 127, 255}));
      connect(tailRadial.port_a, tailValveOut.port_a) annotation(
        Line(points = {{160, 12}, {160, 36}, {184, 36}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[3], bodyValveIn.pin_p) annotation(
        Line(points = {{-198, 130}, {-52, 130}, {-52, 50}, {-46, 50}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[4], bodyValveOut.pin_p) annotation(
        Line(points = {{-198, 130}, {16, 130}, {16, 50}, {28, 50}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[5], tailValveIn.pin_p) annotation(
        Line(points = {{-198, 130}, {108, 130}, {108, 50}, {116, 50}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[6], tailValveOut.pin_p) annotation(
        Line(points = {{-198, 130}, {178, 130}, {178, 50}, {184, 50}}, color = {255, 0, 0}));
      connect(bodyValveIn.pin_n, ground[3].p) annotation(
        Line(points = {{-18, 50}, {-12, 50}, {-12, 100}, {-260, 100}, {-260, 115}}, color = {255, 0, 0}));
      connect(bodyValveOut.pin_n, ground[4].p) annotation(
        Line(points = {{56, 50}, {70, 50}, {70, 100}, {-260, 100}, {-260, 115}}, color = {255, 0, 0}));
      connect(tailValveIn.pin_n, ground[5].p) annotation(
        Line(points = {{144, 50}, {152, 50}, {152, 100}, {-260, 100}, {-260, 115}}, color = {255, 0, 0}));
      connect(tailValveOut.pin_n, ground[6].p) annotation(
        Line(points = {{212, 50}, {218, 50}, {218, 100}, {-260, 100}, {-260, 115}}, color = {255, 0, 0}));
      connect(headValveOut.port_b, ambient.ports[1]) annotation(
        Line(points = {{-82, 36}, {-68, 36}, {-68, 63}, {232, 63}}, color = {0, 127, 255}));
      connect(bodyValveOut.port_b, ambient.ports[2]) annotation(
        Line(points = {{56, 36}, {84, 36}, {84, 63}, {232, 63}}, color = {0, 127, 255}));
      connect(tailValveOut.port_b, ambient.ports[3]) annotation(
        Line(points = {{212, 36}, {222, 36}, {222, 63}, {232, 63}}, color = {0, 127, 255}));
      connect(headValveOut.port_a, headRadial.port_a) annotation(
        Line(points = {{-110, 36}, {-140, 36}, {-140, 12}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[1], headValveIn.pin_p) annotation(
        Line(points = {{-198, 130}, {-186, 130}, {-186, 50}, {-176, 50}}, color = {255, 0, 0}));
      airIn.onOff[1] = lUTControl.controllerLUT.y[1];
      airIn.onOff[2] = lUTControl.controllerLUT.y[3];
      airIn.onOff[3] = lUTControl.controllerLUT.y[5];
      connect(tailRadial.flange_a, bodyLinear.flange_b) annotation(
        Line(points = {{140, -8}, {20, -8}}, color = {0, 127, 0}));
      connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
        Line(points = {{-20, -8}, {-120, -8}}, color = {0, 127, 0}));
      connect(airIn.port_b[1], headValveIn.port_a) annotation(
        Line(points = {{-222, 72}, {-200, 72}, {-200, 36}, {-176, 36}}, color = {0, 127, 255}));
      connect(airIn.port_b[2], bodyValveIn.port_a) annotation(
        Line(points = {{-222, 72}, {-60, 72}, {-60, 36}, {-46, 36}}, color = {0, 127, 255}));
      connect(airIn.port_b[3], tailValveIn.port_a) annotation(
        Line(points = {{-222, 72}, {100, 72}, {100, 36}, {116, 36}}, color = {0, 127, 255}));
      connect(headRadial.flange_a, headMass.flange_a) annotation(
        Line(points = {{-160, -8}, {-176, -8}}, color = {0, 127, 0}));
  connect(tailMass.flange_b, tailRadial.flange_b) annotation(
        Line(points = {{200, -8}, {180, -8}}, color = {0, 127, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 75, Tolerance = 1e-06, Interval = 0.001),
        Diagram(coordinateSystem(extent = {{-300, 160}, {280, -40}})),
        Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
    end WormDynamicComp;
    
    model WormFinal
      replaceable package Medium = Modelica.Media.Air.DryAirNasa;
      //
      // Parameters
      //
      parameter Modelica.Units.SI.Mass m_head = 0.4782;
      parameter Modelica.Units.SI.Mass m_tail = 0.4942;
      parameter Modelica.Units.SI.Mass m_l = 0.1681;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 100000 "Nominal pressure difference across the valve";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.0015 "Nominal mass flow rate across the valve at dp_nominal"; // 0.0015
      parameter Modelica.Units.SI.MassFlowRate m_flow_regulator = 0.0000175 "Mass flow through the regulator";       // 0.000015 gives 32cm in 5 cycles
      //
      // Components
      //
      Modelica.Mechanics.Translational.Components.Mass mass(m = m_head) annotation(
        Placement(visible = true, transformation(origin = {-196, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, p_ambient = 100567) annotation(
        Placement(visible = true, transformation(origin = {-70, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.CylindricalActuators.LinearActuatorStaticLS bodyLinear(redeclare package Medium = Medium, l_0 = 0.146, p_ref = system.p_ambient, r_0 = 0.0375) annotation(
        Placement(visible = true, transformation(origin = {1.77636e-15, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.RadialActuators.RadialPressureLS tailRadial(redeclare package Medium = Medium, m = m_tail, p_start = system.p_ambient, r_joint = 0.095) annotation(
        Placement(visible = true, transformation(origin = {160, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3, p = system.p_ambient) annotation(
        Placement(visible = true, transformation(origin = {250, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      WormBot.Components.RadialActuators.RadialPressureLS headRadial(redeclare package Medium = Medium, m = m_head, p_ambient = system.p_ambient, r_joint = 0.095) annotation(
        Placement(visible = true, transformation(origin = {-140, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      WormBot.Components.ValveControl.LUTControl lUTControl annotation(
        Placement(visible = true, transformation(origin = {-210, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-270, 110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      WormBot.Components.Valves.SolenoidValve headValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-166, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve headValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-98, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {36, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve bodyValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {198, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.Valves.SolenoidValve tailValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
        Placement(visible = true, transformation(origin = {130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      WormBot.Components.AirIn airIn(m_flow_regulator = m_flow_regulator) annotation(
        Placement(visible = true, transformation(origin = {-250, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Translational.Components.Mass mass1(m = m_tail) annotation(
        Placement(visible = true, transformation(origin = {212, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
        Line(points = {{-20, -8}, {-120, -8}}, color = {0, 127, 0}));
      connect(headRadial.flange_a, mass.flange_a) annotation(
        Line(points = {{-160, -8}, {-186, -8}}, color = {0, 127, 0}));
      connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
        Line(points = {{20, -8}, {140, -8}}, color = {0, 127, 0}));
      connect(lUTControl.pin_n, ground.p) annotation(
        Line(points = {{-220, 120.2}, {-260, 120.2}, {-260, 110.2}}, color = {255, 0, 0}, thickness = 0.5));
      connect(headValveIn.port_b, headRadial.port_a) annotation(
        Line(points = {{-156, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(headValveIn.pin_n, ground[1].p) annotation(
        Line(points = {{-156, 44}, {-136, 44}, {-136, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[2], headValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {-116, 120}, {-116, 44}, {-108, 44}}, color = {255, 0, 0}));
      connect(headValveOut.pin_n, ground[2].p) annotation(
        Line(points = {{-88, 44}, {-76, 44}, {-76, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(bodyValveIn.port_b, bodyLinear.port_a) annotation(
        Line(points = {{-22, 34}, {0, 34}, {0, 12}}, color = {0, 127, 255}));
      connect(bodyLinear.port_a, bodyValveOut.port_a) annotation(
        Line(points = {{1.77636e-15, 12}, {1.77636e-15, 34}, {26, 34}}, color = {0, 127, 255}));
      connect(tailValveIn.port_b, tailRadial.port_a) annotation(
        Line(points = {{140, 34}, {160, 34}, {160, 12}}, color = {0, 127, 255}));
      connect(tailRadial.port_a, tailValveOut.port_a) annotation(
        Line(points = {{160, 12}, {160, 34}, {188, 34}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[3], bodyValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-52, 120}, {-52, 44}, {-42, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[4], bodyValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {16, 120}, {16, 44}, {26, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[5], tailValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {108, 120}, {108, 44}, {120, 44}}, color = {255, 0, 0}));
      connect(lUTControl.pin_p[6], tailValveOut.pin_p) annotation(
        Line(points = {{-200, 120}, {178, 120}, {178, 44}, {188, 44}}, color = {255, 0, 0}));
      connect(bodyValveIn.pin_n, ground[3].p) annotation(
        Line(points = {{-22, 44}, {-12, 44}, {-12, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(bodyValveOut.pin_n, ground[4].p) annotation(
        Line(points = {{46, 44}, {70, 44}, {70, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveIn.pin_n, ground[5].p) annotation(
        Line(points = {{140, 44}, {152, 44}, {152, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(tailValveOut.pin_n, ground[6].p) annotation(
        Line(points = {{208, 44}, {216, 44}, {216, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
      connect(headValveOut.port_b, ambient.ports[1]) annotation(
        Line(points = {{-88, 34}, {-68, 34}, {-68, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(bodyValveOut.port_b, ambient.ports[2]) annotation(
        Line(points = {{46, 34}, {84, 34}, {84, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(tailValveOut.port_b, ambient.ports[3]) annotation(
        Line(points = {{208, 34}, {222, 34}, {222, 60}, {240, 60}}, color = {0, 127, 255}));
      connect(headValveOut.port_a, headRadial.port_a) annotation(
        Line(points = {{-108, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
      connect(lUTControl.pin_p[1], headValveIn.pin_p) annotation(
        Line(points = {{-200, 120}, {-186, 120}, {-186, 44}, {-176, 44}}, color = {255, 0, 0}));
      airIn.onOff[1] = lUTControl.controllerLUT.y[1];
      airIn.onOff[2] = lUTControl.controllerLUT.y[3];
      airIn.onOff[3] = lUTControl.controllerLUT.y[5];
      connect(airIn.port_b[1], headValveIn.port_a) annotation(
        Line(points = {{-240, 70}, {-200, 70}, {-200, 34}, {-176, 34}}, color = {0, 127, 255}));
      connect(airIn.port_b[2], bodyValveIn.port_a) annotation(
        Line(points = {{-240, 70}, {-60, 70}, {-60, 34}, {-42, 34}}, color = {0, 127, 255}));
      connect(airIn.port_b[3], tailValveIn.port_a) annotation(
        Line(points = {{-240, 70}, {100, 70}, {100, 34}, {120, 34}}, color = {0, 127, 255}));
    connect(mass1.flange_b, tailRadial.flange_b) annotation(
        Line(points = {{202, -8}, {180, -8}}, color = {0, 127, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 75, Tolerance = 1e-06, Interval = 0.01),
        Diagram(coordinateSystem(extent = {{-280, 200}, {400, -40}})),
        Documentation(info = "<html><head></head><body>This is the final complete worm model. This model reads in the .txt file for the worm control, so the filepath to that file may need to be adjusted on difference devices. The linear actuator uses the least-squares polynomial fit model for the elongation</body></html>"));
    end WormFinal;
  end Systems;

  package Tests
    package MechanicalTests
      model FrictionAndOnOffTest
        parameter Modelica.Units.SI.Mass m = 1;
        Modelica.Mechanics.Translational.Sources.Force force annotation(
          Placement(visible = true, transformation(origin = {60, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp(duration = 0.5, height = 3, offset = 0, startTime = 1) annotation(
          Placement(visible = true, transformation(origin = {70, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.Sticking frictionActuator annotation(
          Placement(visible = true, transformation(origin = {20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
          Placement(visible = true, transformation(origin = {-50, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression realExpression(y = m*9.81) annotation(
          Placement(visible = true, transformation(origin = {-80, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.Sticking frictionJoint(mu_k = 0.2, mu_s = 0.25) annotation(
          Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.BooleanStep booleanStep(startTime = 2) annotation(
          Placement(visible = true, transformation(origin = {-80, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.ContactForces.OnOffContact onOffContact annotation(
          Placement(visible = true, transformation(origin = {-30, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(frictionActuator.flange_b, force.flange) annotation(
          Line(points = {{30, 0}, {50, 0}}, color = {0, 127, 0}));
        connect(frictionActuator.flange_a, frictionJoint.flange_b) annotation(
          Line(points = {{10, 0}, {-10, 0}}, color = {0, 127, 0}));
        connect(realExpression.y, onOffContact.fn) annotation(
          Line(points = {{-69, 80}, {-59, 80}, {-59, 56}, {-41, 56}}, color = {0, 0, 127}));
        connect(booleanStep.y, onOffContact.switching) annotation(
          Line(points = {{-69, 40}, {-69, 42}, {-61, 42}, {-61, 44}, {-41, 44}}, color = {255, 0, 255}));
        connect(onOffContact.actuator, frictionActuator.fn) annotation(
          Line(points = {{-20, 56}, {20, 56}, {20, 10}}, color = {0, 0, 127}));
        connect(onOffContact.joint, frictionJoint.fn) annotation(
          Line(points = {{-20, 44}, {0, 44}, {0, 20}, {-20, 20}, {-20, 10}}, color = {0, 0, 127}));
        connect(ramp.y, force.f) annotation(
          Line(points = {{81, 30}, {99, 30}, {99, 0}, {71, 0}}, color = {0, 0, 127}));
        connect(frictionJoint.flange_a, mass.flange_a) annotation(
          Line(points = {{-30, 0}, {-40, 0}}, color = {0, 127, 0}));
        annotation(
          experiment(StartTime = 0, StopTime = 3, Tolerance = 1e-06, Interval = 0.006));
      end FrictionAndOnOffTest;

      model LinearPressureTest
        replaceable package Medium = Modelica.Media.Air.DryAirNasa;
        Components.CylindricalActuators.LinearPressure linearPressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Components.Mass mass(m = 1, v(fixed = true, start = 0)) annotation(
          Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
          Placement(visible = true, transformation(origin = {70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Components.Fixed fixed annotation(
          Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 1, p = Medium.p_default*1.5) annotation(
          Placement(visible = true, transformation(origin = {-90, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {-40, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp(duration = 0.5, offset = 0, startTime = 0.25) annotation(
          Placement(visible = true, transformation(origin = {-90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(mass.flange_b, linearPressure.flange_a) annotation(
          Line(points = {{-40, 0}, {-10, 0}}, color = {0, 127, 0}));
        connect(fixed.flange, linearPressure.flange_b) annotation(
          Line(points = {{40, 0}, {10, 0}}, color = {0, 127, 0}));
        connect(pressurizedAir.ports[1], valveLinear.port_a) annotation(
          Line(points = {{-80, 50}, {-50, 50}}, color = {0, 127, 255}));
        connect(ramp.y, valveLinear.opening) annotation(
          Line(points = {{-78, 90}, {-40, 90}, {-40, 58}}, color = {0, 0, 127}));
        connect(valveLinear.port_b, linearPressure.port_a) annotation(
          Line(points = {{-30, 50}, {0, 50}, {0, 10}}, color = {0, 127, 255}));
      end LinearPressureTest;

      model LinearAndFrictionOnOff
        replaceable package Medium = Modelica.Media.Air.DryAirNasa;
        parameter Modelica.Units.SI.Mass m = 5;
        WormBot.Components.ContactForces.OnOffContact onOffContact annotation(
          Placement(visible = true, transformation(origin = {-50, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
          Placement(visible = true, transformation(origin = {-90, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        WormBot.Components.Friction.LuGre frictionActuator annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.BooleanStep control(startTime = 1) annotation(
          Placement(visible = true, transformation(origin = {-100, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.LuGre frictionJoint(mu_k = 0.2, mu_s = 0.25) annotation(
          Placement(visible = true, transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression normalForce(y = m*9.81) annotation(
          Placement(visible = true, transformation(origin = {-100, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
          Placement(visible = true, transformation(origin = {-30, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 1, p = Medium.p_default*1.5) annotation(
          Placement(visible = true, transformation(origin = {10, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.CylindricalActuators.LinearPressure linearPressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp(duration = 0.5, offset = 0, startTime = 1) annotation(
          Placement(visible = true, transformation(origin = {-10, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Components.Fixed fixed annotation(
          Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear1(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 1) annotation(
          Placement(visible = true, transformation(origin = {130, 70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp1(duration = 0.5, height = 1, offset = 0, startTime = 3) annotation(
          Placement(visible = true, transformation(origin = {70, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp2(duration = 0.5, height = -1, offset = 0, startTime = 2) annotation(
          Placement(visible = true, transformation(origin = {-10, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Math.Add add annotation(
          Placement(visible = true, transformation(origin = {30, 124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(normalForce.y, onOffContact.fn) annotation(
          Line(points = {{-89, 80}, {-79, 80}, {-79, 56}, {-61, 56}}, color = {0, 0, 127}));
        connect(onOffContact.actuator, frictionActuator.fn) annotation(
          Line(points = {{-40, 56}, {0, 56}, {0, 10}}, color = {0, 0, 127}));
        connect(control.y, onOffContact.switching) annotation(
          Line(points = {{-89, 40}, {-89, 42}, {-81, 42}, {-81, 44}, {-61, 44}}, color = {255, 0, 255}));
        connect(onOffContact.joint, frictionJoint.fn) annotation(
          Line(points = {{-40, 44}, {-20, 44}, {-20, 20}, {-40, 20}, {-40, 10}}, color = {0, 0, 127}));
        connect(frictionActuator.flange_a, frictionJoint.flange_b) annotation(
          Line(points = {{-10, 0}, {-30, 0}}, color = {0, 127, 0}));
        connect(pressurizedAir.ports[1], valveLinear.port_a) annotation(
          Line(points = {{20, 70}, {38, 70}}, color = {0, 127, 255}));
        connect(fixed.flange, linearPressure.flange_b) annotation(
          Line(points = {{100, 0}, {80, 0}}, color = {0, 127, 0}));
        connect(valveLinear.port_b, linearPressure.port_a) annotation(
          Line(points = {{58, 70}, {70, 70}, {70, 10}}, color = {0, 127, 255}));
        connect(linearPressure.flange_a, frictionActuator.flange_b) annotation(
          Line(points = {{60, 0}, {10, 0}}, color = {0, 127, 0}));
        connect(frictionJoint.flange_a, mass.flange_a) annotation(
          Line(points = {{-50, 0}, {-80, 0}}, color = {0, 127, 0}));
        connect(valveLinear1.port_a, linearPressure.port_a) annotation(
          Line(points = {{80, 70}, {70, 70}, {70, 10}}, color = {0, 127, 255}));
        connect(valveLinear1.port_b, ambient.ports[1]) annotation(
          Line(points = {{100, 70}, {120, 70}}, color = {0, 127, 255}));
        connect(ramp1.y, valveLinear1.opening) annotation(
          Line(points = {{82, 130}, {90, 130}, {90, 78}}, color = {0, 0, 127}));
        connect(ramp.y, add.u2) annotation(
          Line(points = {{2, 110}, {18, 110}, {18, 118}}, color = {0, 0, 127}));
        connect(ramp2.y, add.u1) annotation(
          Line(points = {{2, 140}, {18, 140}, {18, 130}}, color = {0, 0, 127}));
        connect(add.y, valveLinear.opening) annotation(
          Line(points = {{42, 124}, {48, 124}, {48, 78}}, color = {0, 0, 127}));
      end LinearAndFrictionOnOff;

      model FrictionTest
        parameter Modelica.Units.SI.Mass m = 1;
        Modelica.Mechanics.Translational.Sources.Force force annotation(
          Placement(visible = true, transformation(origin = {40, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp(duration = 1, height = 3, offset = 0, startTime = 0) annotation(
          Placement(visible = true, transformation(origin = {50, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
          Placement(visible = true, transformation(origin = {-38, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression realExpression(y = m*9.81) annotation(
          Placement(visible = true, transformation(origin = {-28, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.Sticking frictionJoint(mu_k = 0.2, mu_s = 0.25) annotation(
          Placement(visible = true, transformation(origin = {-8, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(ramp.y, force.f) annotation(
          Line(points = {{61, 60}, {79, 60}, {79, 30}, {51, 30}}, color = {0, 0, 127}));
        connect(frictionJoint.flange_a, mass.flange_a) annotation(
          Line(points = {{-18, 30}, {-28, 30}}, color = {0, 127, 0}));
        connect(realExpression.y, frictionJoint.fn) annotation(
          Line(points = {{-17, 54}, {-8, 54}, {-8, 40}}, color = {0, 0, 127}));
        connect(force.flange, frictionJoint.flange_b) annotation(
          Line(points = {{30, 30}, {2, 30}}, color = {0, 127, 0}));
        annotation(
          experiment(StartTime = 0, StopTime = 3, Tolerance = 1e-06, Interval = 0.006));
      end FrictionTest;
    end MechanicalTests;

    package MechAndFluidTests
      model RadialAndLinearTest
        replaceable package Medium = Modelica.Media.Air.DryAirNasa;
        parameter Modelica.Units.SI.Mass m = 5;
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
          Placement(visible = true, transformation(origin = {-90, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
          Placement(visible = true, transformation(origin = {-30, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {48, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 2, p = Medium.p_default*1.5) annotation(
          Placement(visible = true, transformation(origin = {-108, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.CylindricalActuators.LinearPressure linearPressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp(duration = 0.5, offset = 0, startTime = 1) annotation(
          Placement(visible = true, transformation(origin = {-10, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Components.Fixed fixed annotation(
          Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear1(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 1) annotation(
          Placement(visible = true, transformation(origin = {130, 70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp1(duration = 0.5, height = 1, offset = 0, startTime = 3) annotation(
          Placement(visible = true, transformation(origin = {70, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp2(duration = 0.5, height = -1, offset = 0, startTime = 2) annotation(
          Placement(visible = true, transformation(origin = {-10, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Math.Add add annotation(
          Placement(visible = true, transformation(origin = {30, 124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.SphericalActuators.RadialPressure radialPressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {-30, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear2(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp3(duration = 0.5, offset = 0, startTime = 1) annotation(
          Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(pressurizedAir.ports[1], valveLinear.port_a) annotation(
          Line(points = {{-98, 70}, {38, 70}}, color = {0, 127, 255}));
        connect(fixed.flange, linearPressure.flange_b) annotation(
          Line(points = {{100, 0}, {80, 0}}, color = {0, 127, 0}));
        connect(valveLinear.port_b, linearPressure.port_a) annotation(
          Line(points = {{58, 70}, {70, 70}, {70, 10}}, color = {0, 127, 255}));
        connect(valveLinear1.port_a, linearPressure.port_a) annotation(
          Line(points = {{80, 70}, {70, 70}, {70, 10}}, color = {0, 127, 255}));
        connect(valveLinear1.port_b, ambient.ports[1]) annotation(
          Line(points = {{100, 70}, {120, 70}}, color = {0, 127, 255}));
        connect(ramp1.y, valveLinear1.opening) annotation(
          Line(points = {{82, 130}, {90, 130}, {90, 78}}, color = {0, 0, 127}));
        connect(ramp.y, add.u2) annotation(
          Line(points = {{2, 110}, {18, 110}, {18, 118}}, color = {0, 0, 127}));
        connect(ramp2.y, add.u1) annotation(
          Line(points = {{2, 140}, {18, 140}, {18, 130}}, color = {0, 0, 127}));
        connect(add.y, valveLinear.opening) annotation(
          Line(points = {{42, 124}, {48, 124}, {48, 78}}, color = {0, 0, 127}));
        connect(linearPressure.flange_a, radialPressure.flange_b) annotation(
          Line(points = {{60, 0}, {-20, 0}}, color = {0, 127, 0}));
        connect(radialPressure.flange_a, mass.flange_a) annotation(
          Line(points = {{-40, 0}, {-80, 0}}, color = {0, 127, 0}));
        connect(pressurizedAir.ports[2], valveLinear2.port_a) annotation(
          Line(points = {{-98, 70}, {-90, 70}, {-90, 30}, {-60, 30}}, color = {0, 127, 255}));
        connect(valveLinear2.port_b, radialPressure.port_a) annotation(
          Line(points = {{-40, 30}, {-30, 30}, {-30, 10}}, color = {0, 127, 255}));
        connect(ramp3.y, valveLinear2.opening) annotation(
          Line(points = {{-59, 50}, {-50, 50}, {-50, 38}}, color = {0, 0, 127}));
        annotation(
          experiment(StartTime = 0, StopTime = 6, Tolerance = 1e-6, Interval = 0.0001));
      end RadialAndLinearTest;

      model WormTest
        replaceable package Medium = Modelica.Media.Air.DryAirNasa;
        parameter Modelica.Units.SI.Mass m = 1;
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
          Placement(visible = true, transformation(origin = {-120, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
          Placement(visible = true, transformation(origin = {90, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {-22, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 2, p = Medium.p_default*1.5) annotation(
          Placement(visible = true, transformation(origin = {-150, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.CylindricalActuators.LinearPressure bodyLinear(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp(duration = 0.5, offset = 0, startTime = 1) annotation(
          Placement(visible = true, transformation(origin = {-80, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure tailRadial(redeclare package Medium = Medium, p_start = 1.5*system.p_start) annotation(
          Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear1(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {20, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 2) annotation(
          Placement(visible = true, transformation(origin = {130, 70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp1(duration = 0.5, height = 1, offset = 0, startTime = 4) annotation(
          Placement(visible = true, transformation(origin = {0, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp2(duration = 0.5, height = -1, offset = 0, startTime = 2) annotation(
          Placement(visible = true, transformation(origin = {-80, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Math.Add add annotation(
          Placement(visible = true, transformation(origin = {-40, 124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure headRadial(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {-60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear2(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {-80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp3(duration = 0.5, offset = 0, startTime = 3) annotation(
          Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear3(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {86, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp4(duration = 0.5, offset = 0, startTime = 3) annotation(
          Placement(visible = true, transformation(origin = {46, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(pressurizedAir.ports[1], valveLinear.port_a) annotation(
          Line(points = {{-140, 80}, {-32, 80}}, color = {0, 127, 255}));
        connect(valveLinear1.port_b, ambient.ports[1]) annotation(
          Line(points = {{30, 80}, {75, 80}, {75, 70}, {120, 70}}, color = {0, 127, 255}));
        connect(ramp1.y, valveLinear1.opening) annotation(
          Line(points = {{11, 130}, {20, 130}, {20, 88}}, color = {0, 0, 127}));
        connect(ramp.y, add.u2) annotation(
          Line(points = {{-69, 110}, {-53, 110}, {-53, 118}}, color = {0, 0, 127}));
        connect(ramp2.y, add.u1) annotation(
          Line(points = {{-69, 140}, {-53, 140}, {-53, 130}}, color = {0, 0, 127}));
        connect(add.y, valveLinear.opening) annotation(
          Line(points = {{-29, 124}, {-22, 124}, {-22, 88}}, color = {0, 0, 127}));
        connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
          Line(points = {{-10, 0}, {-50, 0}}, color = {0, 127, 0}));
        connect(headRadial.flange_a, mass.flange_a) annotation(
          Line(points = {{-70, 0}, {-110, 0}}, color = {0, 127, 0}));
        connect(valveLinear2.port_b, headRadial.port_a) annotation(
          Line(points = {{-70, 30}, {-60, 30}, {-60, 10}}, color = {0, 127, 255}));
        connect(ramp3.y, valveLinear2.opening) annotation(
          Line(points = {{-89, 50}, {-80, 50}, {-80, 38}}, color = {0, 0, 127}));
        connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
          Line(points = {{10, 0}, {50, 0}}, color = {0, 127, 0}));
        connect(valveLinear.port_b, bodyLinear.port_a) annotation(
          Line(points = {{-12, 80}, {0, 80}, {0, 10}}, color = {0, 127, 255}));
        connect(bodyLinear.port_a, valveLinear1.port_a) annotation(
          Line(points = {{0, 10}, {0, 80}, {10, 80}}, color = {0, 127, 255}));
        connect(ramp4.y, valveLinear3.opening) annotation(
          Line(points = {{57, 58}, {85.5, 58}, {85.5, 48}, {86, 48}}, color = {0, 0, 127}));
        connect(tailRadial.port_a, valveLinear3.port_a) annotation(
          Line(points = {{60, 10}, {60, 40}, {76, 40}}, color = {0, 127, 255}));
        connect(valveLinear3.port_b, ambient.ports[2]) annotation(
          Line(points = {{96, 40}, {108, 40}, {108, 70}, {120, 70}}, color = {0, 127, 255}));
        connect(pressurizedAir.ports[2], valveLinear2.port_a) annotation(
          Line(points = {{-140, 80}, {-120, 80}, {-120, 30}, {-90, 30}}, color = {0, 127, 255}));
        annotation(
          experiment(StartTime = 0, StopTime = 6, Tolerance = 1e-06, Interval = 0.0001),
          Diagram(coordinateSystem(extent = {{-160, 160}, {140, -20}})),
          Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
      end WormTest;

      model WormTestStatic
        replaceable package Medium = Modelica.Media.Air.DryAirNasa;
        parameter Modelica.Units.SI.Mass m = 1;
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
          Placement(visible = true, transformation(origin = {-120, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
          Placement(visible = true, transformation(origin = {90, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {-22, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 2, p = Medium.p_default*2) annotation(
          Placement(visible = true, transformation(origin = {-150, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.CylindricalActuators.LinearActuatorStatic bodyLinear(redeclare package Medium = Medium, L = 0.150, R = 0.0375, t = 0.002) annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp(duration = 0.5, offset = 0, startTime = 1) annotation(
          Placement(visible = true, transformation(origin = {-80, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure tailRadial(redeclare package Medium = Medium, gamma_max = 2, p_start = 2*system.p_start) annotation(
          Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear1(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {20, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 2) annotation(
          Placement(visible = true, transformation(origin = {130, 70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp1(duration = 0.5, height = 1, offset = 0, startTime = 4) annotation(
          Placement(visible = true, transformation(origin = {0, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp2(duration = 0.5, height = -1, offset = 0, startTime = 2) annotation(
          Placement(visible = true, transformation(origin = {-80, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Math.Add add annotation(
          Placement(visible = true, transformation(origin = {-40, 124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure headRadial(redeclare package Medium = Medium, gamma_max = 2) annotation(
          Placement(visible = true, transformation(origin = {-60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear2(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {-80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp3(duration = 0.5, offset = 0, startTime = 3) annotation(
          Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear3(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {86, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp4(duration = 0.5, offset = 0, startTime = 3) annotation(
          Placement(visible = true, transformation(origin = {46, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(pressurizedAir.ports[1], valveLinear.port_a) annotation(
          Line(points = {{-140, 80}, {-32, 80}}, color = {0, 127, 255}));
        connect(valveLinear1.port_b, ambient.ports[1]) annotation(
          Line(points = {{30, 80}, {75, 80}, {75, 70}, {120, 70}}, color = {0, 127, 255}));
        connect(ramp1.y, valveLinear1.opening) annotation(
          Line(points = {{11, 130}, {20, 130}, {20, 88}}, color = {0, 0, 127}));
        connect(ramp.y, add.u2) annotation(
          Line(points = {{-69, 110}, {-53, 110}, {-53, 118}}, color = {0, 0, 127}));
        connect(ramp2.y, add.u1) annotation(
          Line(points = {{-69, 140}, {-53, 140}, {-53, 130}}, color = {0, 0, 127}));
        connect(add.y, valveLinear.opening) annotation(
          Line(points = {{-29, 124}, {-22, 124}, {-22, 88}}, color = {0, 0, 127}));
        connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
          Line(points = {{-10, 0}, {-50, 0}}, color = {0, 127, 0}));
        connect(headRadial.flange_a, mass.flange_a) annotation(
          Line(points = {{-70, 0}, {-110, 0}}, color = {0, 127, 0}));
        connect(valveLinear2.port_b, headRadial.port_a) annotation(
          Line(points = {{-70, 30}, {-60, 30}, {-60, 10}}, color = {0, 127, 255}));
        connect(ramp3.y, valveLinear2.opening) annotation(
          Line(points = {{-89, 50}, {-80, 50}, {-80, 38}}, color = {0, 0, 127}));
        connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
          Line(points = {{10, 0}, {50, 0}}, color = {0, 127, 0}));
        connect(valveLinear.port_b, bodyLinear.port_a) annotation(
          Line(points = {{-12, 80}, {0, 80}, {0, 10}}, color = {0, 127, 255}));
        connect(bodyLinear.port_a, valveLinear1.port_a) annotation(
          Line(points = {{0, 10}, {0, 80}, {10, 80}}, color = {0, 127, 255}));
        connect(ramp4.y, valveLinear3.opening) annotation(
          Line(points = {{57, 58}, {85.5, 58}, {85.5, 48}, {86, 48}}, color = {0, 0, 127}));
        connect(tailRadial.port_a, valveLinear3.port_a) annotation(
          Line(points = {{60, 10}, {60, 40}, {76, 40}}, color = {0, 127, 255}));
        connect(valveLinear3.port_b, ambient.ports[2]) annotation(
          Line(points = {{96, 40}, {108, 40}, {108, 70}, {120, 70}}, color = {0, 127, 255}));
        connect(pressurizedAir.ports[2], valveLinear2.port_a) annotation(
          Line(points = {{-140, 80}, {-120, 80}, {-120, 30}, {-90, 30}}, color = {0, 127, 255}));
        annotation(
          experiment(StartTime = 0, StopTime = 6, Tolerance = 1e-06, Interval = 0.0001),
          Diagram(coordinateSystem(extent = {{-160, 160}, {140, -20}})),
          Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
      end WormTestStatic;

      model WormTestDynamic
        replaceable package Medium = Modelica.Media.Air.DryAirNasa;
        parameter Modelica.Units.SI.Mass m = 1;
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
          Placement(visible = true, transformation(origin = {-120, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
          Placement(visible = true, transformation(origin = {90, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {-22, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 2, p = Medium.p_default*2) annotation(
          Placement(visible = true, transformation(origin = {-150, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.CylindricalActuators.LinearActuatorDynamic bodyLinear(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp(duration = 0.5, offset = 0, startTime = 1) annotation(
          Placement(visible = true, transformation(origin = {-80, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure tailRadial(redeclare package Medium = Medium, gamma_max = 2, mu_k_actuator = 1, mu_s_actuator = 1, p_start = 2*system.p_start) annotation(
          Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear1(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {20, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 2) annotation(
          Placement(visible = true, transformation(origin = {130, 70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp1(duration = 0.5, height = 1, offset = 0, startTime = 4) annotation(
          Placement(visible = true, transformation(origin = {0, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp2(duration = 0.5, height = -1, offset = 0, startTime = 2) annotation(
          Placement(visible = true, transformation(origin = {-80, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Math.Add add annotation(
          Placement(visible = true, transformation(origin = {-40, 124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure headRadial(redeclare package Medium = Medium, gamma_max = 2, mu_k_actuator = 1, mu_s_actuator = 1) annotation(
          Placement(visible = true, transformation(origin = {-60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear2(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {-80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp3(duration = 0.5, offset = 0, startTime = 3) annotation(
          Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear3(redeclare package Medium = Medium, dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {86, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp4(duration = 0.5, offset = 0, startTime = 3) annotation(
          Placement(visible = true, transformation(origin = {46, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(pressurizedAir.ports[1], valveLinear.port_a) annotation(
          Line(points = {{-140, 80}, {-32, 80}}, color = {0, 127, 255}));
        connect(valveLinear1.port_b, ambient.ports[1]) annotation(
          Line(points = {{30, 80}, {75, 80}, {75, 70}, {120, 70}}, color = {0, 127, 255}));
        connect(ramp1.y, valveLinear1.opening) annotation(
          Line(points = {{11, 130}, {20, 130}, {20, 88}}, color = {0, 0, 127}));
        connect(ramp.y, add.u2) annotation(
          Line(points = {{-69, 110}, {-53, 110}, {-53, 118}}, color = {0, 0, 127}));
        connect(ramp2.y, add.u1) annotation(
          Line(points = {{-69, 140}, {-53, 140}, {-53, 130}}, color = {0, 0, 127}));
        connect(add.y, valveLinear.opening) annotation(
          Line(points = {{-29, 124}, {-22, 124}, {-22, 88}}, color = {0, 0, 127}));
        connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
          Line(points = {{-10, 0}, {-50, 0}}, color = {0, 127, 0}));
        connect(headRadial.flange_a, mass.flange_a) annotation(
          Line(points = {{-70, 0}, {-110, 0}}, color = {0, 127, 0}));
        connect(valveLinear2.port_b, headRadial.port_a) annotation(
          Line(points = {{-70, 30}, {-60, 30}, {-60, 10}}, color = {0, 127, 255}));
        connect(ramp3.y, valveLinear2.opening) annotation(
          Line(points = {{-89, 50}, {-80, 50}, {-80, 38}}, color = {0, 0, 127}));
        connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
          Line(points = {{10, 0}, {50, 0}}, color = {0, 127, 0}));
        connect(valveLinear.port_b, bodyLinear.port_a) annotation(
          Line(points = {{-12, 80}, {0, 80}, {0, 10}}, color = {0, 127, 255}));
        connect(bodyLinear.port_a, valveLinear1.port_a) annotation(
          Line(points = {{0, 10}, {0, 80}, {10, 80}}, color = {0, 127, 255}));
        connect(ramp4.y, valveLinear3.opening) annotation(
          Line(points = {{57, 58}, {85.5, 58}, {85.5, 48}, {86, 48}}, color = {0, 0, 127}));
        connect(tailRadial.port_a, valveLinear3.port_a) annotation(
          Line(points = {{60, 10}, {60, 40}, {76, 40}}, color = {0, 127, 255}));
        connect(valveLinear3.port_b, ambient.ports[2]) annotation(
          Line(points = {{96, 40}, {108, 40}, {108, 70}, {120, 70}}, color = {0, 127, 255}));
        connect(pressurizedAir.ports[2], valveLinear2.port_a) annotation(
          Line(points = {{-140, 80}, {-120, 80}, {-120, 30}, {-90, 30}}, color = {0, 127, 255}));
        annotation(
          experiment(StartTime = 0, StopTime = 671, Tolerance = 1e-06, Interval = 0.01),
          Diagram(coordinateSystem(extent = {{-160, 160}, {140, -20}})),
          Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
      end WormTestDynamic;

      model WormTestStaticFull
        replaceable package Medium = Modelica.Media.Air.DryAirNasa;
        parameter Modelica.Units.SI.Mass m = 1;
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
          Placement(visible = true, transformation(origin = {-120, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
          Placement(visible = true, transformation(origin = {110, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 2, p = Medium.p_default*2) annotation(
          Placement(visible = true, transformation(origin = {-150, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.CylindricalActuators.LinearActuatorStatic bodyLinear(redeclare package Medium = Medium, L = 0.150, R = 0.0375, t = 0.002) annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure tailRadial(redeclare package Medium = Medium, gamma_max = 2, p_start = 2*system.p_start) annotation(
          Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 2) annotation(
          Placement(visible = true, transformation(origin = {130, 70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure headRadial(redeclare package Medium = Medium, gamma_max = 2) annotation(
          Placement(visible = true, transformation(origin = {-60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.ValveControl.LUTControl lUTControl annotation(
          Placement(visible = true, transformation(origin = {-70, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
          Placement(visible = true, transformation(origin = {-130, 130}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
        Modelica.Electrical.Analog.Sensors.VoltageSensor[6] voltageSensor annotation(
          Placement(visible = true, transformation(origin = {-70, 176}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear3(dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {86, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear1(dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {20, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear(dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {-22, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear2(dp_nominal = 999.9999999999999, m_flow_nominal = 0.001) annotation(
          Placement(visible = true, transformation(origin = {-80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Electrical.Analog.Basic.Resistor resistor(R = 1) annotation(
          Placement(visible = true, transformation(origin = {-72, 106}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
          Line(points = {{-10, 0}, {-50, 0}}, color = {0, 127, 0}));
        connect(headRadial.flange_a, mass.flange_a) annotation(
          Line(points = {{-70, 0}, {-110, 0}}, color = {0, 127, 0}));
        connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
          Line(points = {{10, 0}, {50, 0}}, color = {0, 127, 0}));
        connect(lUTControl.pin_n, ground.p) annotation(
          Line(points = {{-80, 140}, {-120, 140}, {-120, 130}}, color = {0, 0, 255}, thickness = 0.5));
        connect(voltageSensor.n, lUTControl.pin_n) annotation(
          Line(points = {{-80, 176}, {-100, 176}, {-100, 140}, {-80, 140}}, color = {0, 0, 255}, thickness = 0.5));
        connect(voltageSensor.p, lUTControl.pin_p) annotation(
          Line(points = {{-60, 176}, {-40, 176}, {-40, 140}, {-60, 140}}, color = {0, 0, 255}, thickness = 0.5));
        connect(voltageSensor[6].v, valveLinear3.opening) annotation(
          Line(points = {{-70, 188}, {86, 188}, {86, 48}}, color = {0, 0, 127}));
        connect(tailRadial.port_a, valveLinear3.port_a) annotation(
          Line(points = {{60, 10}, {60, 40}, {76, 40}}, color = {0, 127, 255}));
        connect(valveLinear3.port_b, ambient.ports[1]) annotation(
          Line(points = {{96, 40}, {108, 40}, {108, 70}, {120, 70}}, color = {0, 127, 255}));
        connect(voltageSensor[4].v, valveLinear1.opening) annotation(
          Line(points = {{-70, 188}, {20, 188}, {20, 88}}, color = {0, 0, 127}));
        connect(bodyLinear.port_a, valveLinear1.port_a) annotation(
          Line(points = {{0, 10}, {0, 80}, {10, 80}}, color = {0, 127, 255}));
        connect(valveLinear1.port_b, ambient.ports[2]) annotation(
          Line(points = {{30, 80}, {75, 80}, {75, 70}, {120, 70}}, color = {0, 127, 255}));
        connect(voltageSensor[3].v, valveLinear.opening) annotation(
          Line(points = {{-70, 188}, {-22, 188}, {-22, 88}}, color = {0, 0, 127}));
        connect(pressurizedAir.ports[1], valveLinear.port_a) annotation(
          Line(points = {{-140, 80}, {-32, 80}}, color = {0, 127, 255}));
        connect(valveLinear.port_b, bodyLinear.port_a) annotation(
          Line(points = {{-12, 80}, {0, 80}, {0, 10}}, color = {0, 127, 255}));
        connect(voltageSensor[1].v, valveLinear2.opening) annotation(
          Line(points = {{-70, 188}, {-106, 188}, {-106, 58}, {-80, 58}, {-80, 38}}, color = {0, 0, 127}));
        connect(valveLinear2.port_b, headRadial.port_a) annotation(
          Line(points = {{-70, 30}, {-60, 30}, {-60, 10}}, color = {0, 127, 255}));
        connect(pressurizedAir.ports[2], valveLinear2.port_a) annotation(
          Line(points = {{-140, 80}, {-120, 80}, {-120, 30}, {-90, 30}}, color = {0, 127, 255}));
        connect(resistor.n, lUTControl.pin_p) annotation(
          Line(points = {{-62, 106}, {-40, 106}, {-40, 140}, {-60, 140}}, color = {0, 0, 255}, thickness = 0.5));
        connect(ground.p, resistor.p) annotation(
          Line(points = {{-120, 130}, {-120, 106}, {-82, 106}}, color = {0, 0, 255}, thickness = 0.5));
        annotation(
          experiment(StartTime = 0, StopTime = 6, Tolerance = 1e-06, Interval = 0.0001),
          Diagram(coordinateSystem(extent = {{-160, 160}, {140, -20}})),
          Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
      end WormTestStaticFull;

      model WormTestStaticFullValves
        replaceable package Medium = Modelica.Media.Air.DryAirNasa;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Mass m = 1;
        parameter Modelica.Units.SI.PressureDifference dp_nominal = 100000 "Nominal pressure difference across the valve";
        parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.00237 "Nominal mass flow rate across the valve at dp_nominal";
        //
        // Components
        //
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
          Placement(visible = true, transformation(origin = {-196, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
          Placement(visible = true, transformation(origin = {-70, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 3, p = Medium.p_default*1.5) annotation(
          Placement(visible = true, transformation(origin = {-226, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.CylindricalActuators.LinearActuatorStatic bodyLinear(redeclare package Medium = Medium, L = 0.150, R = 0.0375, t = 0.002) annotation(
          Placement(visible = true, transformation(origin = {1.77636e-15, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure tailRadial(redeclare package Medium = Medium, gamma_max = 2, p_start = 2*system.p_start) annotation(
          Placement(visible = true, transformation(origin = {160, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3) annotation(
          Placement(visible = true, transformation(origin = {236, 62}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure headRadial(redeclare package Medium = Medium, gamma_max = 2) annotation(
          Placement(visible = true, transformation(origin = {-140, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        WormBot.Components.ValveControl.LUTControl lUTControl annotation(
          Placement(visible = true, transformation(origin = {-210, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
          Placement(visible = true, transformation(origin = {-270, 110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
        WormBot.Components.Valves.SolenoidValve headValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {-166, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Valves.SolenoidValve headValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {-98, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Valves.SolenoidValve bodyValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {36, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Valves.SolenoidValve bodyValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Valves.SolenoidValve tailValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {198, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Valves.SolenoidValve tailValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
          Line(points = {{-20, -8}, {-120, -8}}, color = {0, 127, 0}));
        connect(headRadial.flange_a, mass.flange_a) annotation(
          Line(points = {{-160, -8}, {-186, -8}}, color = {0, 127, 0}));
        connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
          Line(points = {{20, -8}, {140, -8}}, color = {0, 127, 0}));
        connect(lUTControl.pin_n, ground.p) annotation(
          Line(points = {{-220, 120.2}, {-260, 120.2}, {-260, 110.2}}, color = {255, 0, 0}, thickness = 0.5));
        connect(headValveIn.port_b, headRadial.port_a) annotation(
          Line(points = {{-156, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
        connect(headValveIn.port_a, pressurizedAir.ports[1]) annotation(
          Line(points = {{-176, 34}, {-196, 34}, {-196, 72}, {-216, 72}}, color = {0, 127, 255}));
        connect(headValveIn.pin_n, ground[1].p) annotation(
          Line(points = {{-156, 44}, {-136, 44}, {-136, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
        connect(lUTControl.pin_p[2], headValveOut.pin_p) annotation(
          Line(points = {{-200, 120}, {-116, 120}, {-116, 44}, {-108, 44}}, color = {255, 0, 0}));
        connect(headValveOut.pin_n, ground[2].p) annotation(
          Line(points = {{-88, 44}, {-76, 44}, {-76, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
        connect(pressurizedAir.ports[2], bodyValveIn.port_a) annotation(
          Line(points = {{-216, 72}, {-60, 72}, {-60, 34}, {-42, 34}}, color = {0, 127, 255}));
        connect(tailValveIn.port_a, pressurizedAir.ports[3]) annotation(
          Line(points = {{120, 34}, {100, 34}, {100, 72}, {-216, 72}}, color = {0, 127, 255}));
        connect(bodyValveIn.port_b, bodyLinear.port_a) annotation(
          Line(points = {{-22, 34}, {0, 34}, {0, 12}}, color = {0, 127, 255}));
        connect(bodyLinear.port_a, bodyValveOut.port_a) annotation(
          Line(points = {{1.77636e-15, 12}, {1.77636e-15, 34}, {26, 34}}, color = {0, 127, 255}));
        connect(tailValveIn.port_b, tailRadial.port_a) annotation(
          Line(points = {{140, 34}, {160, 34}, {160, 12}}, color = {0, 127, 255}));
        connect(tailRadial.port_a, tailValveOut.port_a) annotation(
          Line(points = {{160, 12}, {160, 34}, {188, 34}}, color = {0, 127, 255}));
        connect(lUTControl.pin_p[3], bodyValveIn.pin_p) annotation(
          Line(points = {{-200, 120}, {-52, 120}, {-52, 44}, {-42, 44}}, color = {255, 0, 0}));
        connect(lUTControl.pin_p[4], bodyValveOut.pin_p) annotation(
          Line(points = {{-200, 120}, {16, 120}, {16, 44}, {26, 44}}, color = {255, 0, 0}));
        connect(lUTControl.pin_p[5], tailValveIn.pin_p) annotation(
          Line(points = {{-200, 120}, {108, 120}, {108, 44}, {120, 44}}, color = {255, 0, 0}));
        connect(lUTControl.pin_p[6], tailValveOut.pin_p) annotation(
          Line(points = {{-200, 120}, {178, 120}, {178, 44}, {188, 44}}, color = {255, 0, 0}));
        connect(bodyValveIn.pin_n, ground[3].p) annotation(
          Line(points = {{-22, 44}, {-12, 44}, {-12, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
        connect(bodyValveOut.pin_n, ground[4].p) annotation(
          Line(points = {{46, 44}, {70, 44}, {70, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
        connect(tailValveIn.pin_n, ground[5].p) annotation(
          Line(points = {{140, 44}, {152, 44}, {152, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
        connect(tailValveOut.pin_n, ground[6].p) annotation(
          Line(points = {{208, 44}, {216, 44}, {216, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
        connect(headValveOut.port_b, ambient.ports[1]) annotation(
          Line(points = {{-88, 34}, {-68, 34}, {-68, 62}, {226, 62}}, color = {0, 127, 255}));
        connect(bodyValveOut.port_b, ambient.ports[2]) annotation(
          Line(points = {{46, 34}, {84, 34}, {84, 62}, {226, 62}}, color = {0, 127, 255}));
        connect(tailValveOut.port_b, ambient.ports[3]) annotation(
          Line(points = {{208, 34}, {222, 34}, {222, 62}, {226, 62}}, color = {0, 127, 255}));
        connect(headValveOut.port_a, headRadial.port_a) annotation(
          Line(points = {{-108, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
        connect(lUTControl.pin_p[1], headValveIn.pin_p) annotation(
          Line(points = {{-200, 120}, {-186, 120}, {-186, 44}, {-176, 44}}, color = {255, 0, 0}));
        annotation(
          experiment(StartTime = 0, StopTime = 350, Tolerance = 1e-06, Interval = 0.01),
          Diagram(coordinateSystem(extent = {{-280, 180}, {260, -40}})),
          Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
      end WormTestStaticFullValves;

      model WormTestStaticFullValves1
        replaceable package Medium = Modelica.Media.Air.DryAirNasa;
        parameter Modelica.Units.SI.Mass m = 1;
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
          Placement(visible = true, transformation(origin = {-196, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
          Placement(visible = true, transformation(origin = {76, 156}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 1, p = Medium.p_default*2) annotation(
          Placement(visible = true, transformation(origin = {-226, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.CylindricalActuators.LinearActuatorStatic bodyLinear(redeclare package Medium = Medium, L = 0.150, R = 0.0375, t = 0.002) annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure tailRadial(redeclare package Medium = Medium, gamma_max = 2, p_start = 2*system.p_start) annotation(
          Placement(visible = true, transformation(origin = {166, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3) annotation(
          Placement(visible = true, transformation(origin = {236, 70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure headRadial(redeclare package Medium = Medium, gamma_max = 2) annotation(
          Placement(visible = true, transformation(origin = {-136, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.ValveControl.LUTControl lUTControl annotation(
          Placement(visible = true, transformation(origin = {-146, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
          Placement(visible = true, transformation(origin = {-206, 130}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
        Modelica.Electrical.Analog.Sensors.VoltageSensor[6] voltageSensor annotation(
          Placement(visible = true, transformation(origin = {-146, 176}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
        Components.Valves.SolenoidValve solenoidValve1 annotation(
          Placement(visible = true, transformation(origin = {-98, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Valves.SolenoidValve solenoidValve annotation(
          Placement(visible = true, transformation(origin = {-166, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
          Line(points = {{-10, 0}, {-126, 0}}, color = {0, 127, 0}));
        connect(headRadial.flange_a, mass.flange_a) annotation(
          Line(points = {{-146, 0}, {-186, 0}}, color = {0, 127, 0}));
        connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
          Line(points = {{10, 0}, {156, 0}}, color = {0, 127, 0}));
        connect(lUTControl.pin_n, ground.p) annotation(
          Line(points = {{-156, 140.2}, {-196, 140.2}, {-196, 130.2}}, color = {0, 0, 255}, thickness = 0.5));
        connect(voltageSensor.n, lUTControl.pin_n) annotation(
          Line(points = {{-156, 176}, {-176, 176}, {-176, 140}, {-156, 140}}, color = {0, 0, 255}, thickness = 0.5));
        connect(voltageSensor.p, lUTControl.pin_p) annotation(
          Line(points = {{-136, 176}, {-116, 176}, {-116, 140}, {-136, 140}}, color = {0, 0, 255}, thickness = 0.5));
        connect(tailRadial.port_a, ambient.ports[1]) annotation(
          Line(points = {{166, 10}, {166, 70}, {226, 70}}, color = {0, 127, 255}));
        connect(bodyLinear.port_a, ambient.ports[2]) annotation(
          Line(points = {{0, 10}, {0, 70}, {226, 70}}, color = {0, 127, 255}));
        connect(solenoidValve1.pin_n, ground[2].p) annotation(
          Line(points = {{-88, 52}, {-76, 52}, {-76, 90}, {-196, 90}, {-196, 130}}, color = {0, 0, 255}));
        connect(lUTControl.pin_p[2], solenoidValve1.pin_p) annotation(
          Line(points = {{-136, 140}, {-116, 140}, {-116, 52}, {-108, 52}}, color = {0, 0, 255}));
        connect(headRadial.port_a, solenoidValve1.port_a) annotation(
          Line(points = {{-136, 10}, {-136, 42}, {-108, 42}}, color = {0, 127, 255}));
        connect(solenoidValve.pin_n, ground[1].p) annotation(
          Line(points = {{-156, 52}, {-136, 52}, {-136, 90}, {-196, 90}, {-196, 130}}, color = {0, 0, 255}));
        connect(voltageSensor[1].p, solenoidValve.pin_p) annotation(
          Line(points = {{-136, 176}, {-116, 176}, {-116, 58}, {-180, 58}, {-180, 52}, {-176, 52}}, color = {0, 0, 255}));
        connect(solenoidValve.port_a, pressurizedAir.ports[1]) annotation(
          Line(points = {{-176, 42}, {-196, 42}, {-196, 80}, {-216, 80}}, color = {0, 127, 255}));
        connect(solenoidValve.port_b, headRadial.port_a) annotation(
          Line(points = {{-156, 42}, {-136, 42}, {-136, 10}}, color = {0, 127, 255}));
        connect(solenoidValve1.port_b, ambient.ports[3]) annotation(
          Line(points = {{-88, 42}, {0, 42}, {0, 70}, {226, 70}}, color = {0, 127, 255}));
        annotation(
          experiment(StartTime = 0, StopTime = 60, Tolerance = 1e-06, Interval = 0.12),
          Diagram(coordinateSystem(extent = {{-160, 160}, {140, -20}})),
          Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
      end WormTestStaticFullValves1;

      model WormTestStaticVarVol
        replaceable package Medium = Modelica.Media.Air.DryAirNasa;
        //
        // Parameters
        //
        parameter Modelica.Units.SI.Mass m = 1;
        parameter Modelica.Units.SI.PressureDifference dp_nominal = 100000 "Nominal pressure difference across the valve";
        parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.00237 "Nominal mass flow rate across the valve at dp_nominal";
        //
        // Components
        //
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
          Placement(visible = true, transformation(origin = {-196, -8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
          Placement(visible = true, transformation(origin = {-70, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary pressurizedAir(redeclare package Medium = Medium, nPorts = 3, p = Medium.p_default*1.5) annotation(
          Placement(visible = true, transformation(origin = {-226, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.CylindricalActuators.LinearActuatorStaticVarVol bodyLinear(redeclare package Medium = Medium, L = 0.150, R = 0.0375, t = 0.002) annotation(
          Placement(visible = true, transformation(origin = {1.77636e-15, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure tailRadial(redeclare package Medium = Medium, gamma_max = 2, p_start = 2*system.p_start) annotation(
          Placement(visible = true, transformation(origin = {160, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 3) annotation(
          Placement(visible = true, transformation(origin = {236, 62}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        WormBot.Components.SphericalActuators.RadialPressure headRadial(redeclare package Medium = Medium, gamma_max = 2) annotation(
          Placement(visible = true, transformation(origin = {-140, -8}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        WormBot.Components.ValveControl.LUTControl lUTControl annotation(
          Placement(visible = true, transformation(origin = {-210, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
          Placement(visible = true, transformation(origin = {-270, 110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
        WormBot.Components.Valves.SolenoidValve headValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {-166, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Valves.SolenoidValve headValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {-98, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Valves.SolenoidValve bodyValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {36, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Valves.SolenoidValve bodyValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Valves.SolenoidValve tailValveOut(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {198, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Valves.SolenoidValve tailValveIn(redeclare package Medium = Medium, dp_nominal = dp_nominal, m_flow_nominal = m_flow_nominal) annotation(
          Placement(visible = true, transformation(origin = {130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(bodyLinear.flange_a, headRadial.flange_b) annotation(
          Line(points = {{-20, -8}, {-120, -8}}, color = {0, 127, 0}));
        connect(headRadial.flange_a, mass.flange_a) annotation(
          Line(points = {{-160, -8}, {-186, -8}}, color = {0, 127, 0}));
        connect(bodyLinear.flange_b, tailRadial.flange_a) annotation(
          Line(points = {{20, -8}, {140, -8}}, color = {0, 127, 0}));
        connect(lUTControl.pin_n, ground.p) annotation(
          Line(points = {{-220, 120.2}, {-260, 120.2}, {-260, 110.2}}, color = {255, 0, 0}, thickness = 0.5));
        connect(headValveIn.port_b, headRadial.port_a) annotation(
          Line(points = {{-156, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
        connect(headValveIn.port_a, pressurizedAir.ports[1]) annotation(
          Line(points = {{-176, 34}, {-196, 34}, {-196, 72}, {-216, 72}}, color = {0, 127, 255}));
        connect(headValveIn.pin_n, ground[1].p) annotation(
          Line(points = {{-156, 44}, {-136, 44}, {-136, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
        connect(lUTControl.pin_p[2], headValveOut.pin_p) annotation(
          Line(points = {{-200, 120}, {-116, 120}, {-116, 44}, {-108, 44}}, color = {255, 0, 0}));
        connect(headValveOut.pin_n, ground[2].p) annotation(
          Line(points = {{-88, 44}, {-76, 44}, {-76, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
        connect(pressurizedAir.ports[2], bodyValveIn.port_a) annotation(
          Line(points = {{-216, 72}, {-60, 72}, {-60, 34}, {-42, 34}}, color = {0, 127, 255}));
        connect(tailValveIn.port_a, pressurizedAir.ports[3]) annotation(
          Line(points = {{120, 34}, {100, 34}, {100, 72}, {-216, 72}}, color = {0, 127, 255}));
        connect(bodyValveIn.port_b, bodyLinear.port_a) annotation(
          Line(points = {{-22, 34}, {0, 34}, {0, 12}}, color = {0, 127, 255}));
        connect(bodyLinear.port_a, bodyValveOut.port_a) annotation(
          Line(points = {{1.77636e-15, 12}, {1.77636e-15, 34}, {26, 34}}, color = {0, 127, 255}));
        connect(tailValveIn.port_b, tailRadial.port_a) annotation(
          Line(points = {{140, 34}, {160, 34}, {160, 12}}, color = {0, 127, 255}));
        connect(tailRadial.port_a, tailValveOut.port_a) annotation(
          Line(points = {{160, 12}, {160, 34}, {188, 34}}, color = {0, 127, 255}));
        connect(lUTControl.pin_p[3], bodyValveIn.pin_p) annotation(
          Line(points = {{-200, 120}, {-52, 120}, {-52, 44}, {-42, 44}}, color = {255, 0, 0}));
        connect(lUTControl.pin_p[4], bodyValveOut.pin_p) annotation(
          Line(points = {{-200, 120}, {16, 120}, {16, 44}, {26, 44}}, color = {255, 0, 0}));
        connect(lUTControl.pin_p[5], tailValveIn.pin_p) annotation(
          Line(points = {{-200, 120}, {108, 120}, {108, 44}, {120, 44}}, color = {255, 0, 0}));
        connect(lUTControl.pin_p[6], tailValveOut.pin_p) annotation(
          Line(points = {{-200, 120}, {178, 120}, {178, 44}, {188, 44}}, color = {255, 0, 0}));
        connect(bodyValveIn.pin_n, ground[3].p) annotation(
          Line(points = {{-22, 44}, {-12, 44}, {-12, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
        connect(bodyValveOut.pin_n, ground[4].p) annotation(
          Line(points = {{46, 44}, {70, 44}, {70, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
        connect(tailValveIn.pin_n, ground[5].p) annotation(
          Line(points = {{140, 44}, {152, 44}, {152, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
        connect(tailValveOut.pin_n, ground[6].p) annotation(
          Line(points = {{208, 44}, {216, 44}, {216, 90}, {-260, 90}, {-260, 110}}, color = {255, 0, 0}));
        connect(headValveOut.port_b, ambient.ports[1]) annotation(
          Line(points = {{-88, 34}, {-68, 34}, {-68, 62}, {226, 62}}, color = {0, 127, 255}));
        connect(bodyValveOut.port_b, ambient.ports[2]) annotation(
          Line(points = {{46, 34}, {84, 34}, {84, 62}, {226, 62}}, color = {0, 127, 255}));
        connect(tailValveOut.port_b, ambient.ports[3]) annotation(
          Line(points = {{208, 34}, {222, 34}, {222, 62}, {226, 62}}, color = {0, 127, 255}));
        connect(headValveOut.port_a, headRadial.port_a) annotation(
          Line(points = {{-108, 34}, {-140, 34}, {-140, 12}}, color = {0, 127, 255}));
        connect(lUTControl.pin_p[1], headValveIn.pin_p) annotation(
          Line(points = {{-200, 120}, {-186, 120}, {-186, 44}, {-176, 44}}, color = {255, 0, 0}));
        annotation(
          experiment(StartTime = 0, StopTime = 350, Tolerance = 1e-06, Interval = 0.01),
          Diagram(coordinateSystem(extent = {{-280, 180}, {260, -40}})),
          Documentation(info = "<html><head></head><body>Basic worm-bot simulation. Model completes most of a locomotion cycle. The <b>tailRadial</b>&nbsp;is initialized at full pressure (1.5 atm). with the rest of the components are initialized at 1 atm. At 1 second, <b>valveLinear</b>&nbsp;opens, which pressurizes <b>bodyLinear</b>, moving the head mass forwards. At 2 seconds, <b>valveLinear</b>&nbsp;closes, preventing any more air for entering. At 3 seconds, the <b>tailRadial</b>&nbsp;vents to ambient (depressurizes) and <b>headRadial</b>&nbsp;pressurizes, anchoring the head. At 4 seconds, the <b>bodyLinear</b>&nbsp;is depressurized, pulling to tail.<div><br></div><div>The only part of the locomotion cycle that is not currently modelled is the pressurization of the tail and the depressurization of the head.</div></body></html>"));
      end WormTestStaticVarVol;

      model DynamicTester
        replaceable package Medium = Modelica.Media.Air.DryAirNasa;
        //
        // Components
        //
        Components.CylindricalActuators.LinearActuatorDynamicDirection linearActuatorDynamic1(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {0, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.MassFlowSource_T boundary(redeclare package Medium = Medium, nPorts = 1, use_m_flow_in = true) annotation(
          Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.MassFlowSource_T boundary1(redeclare package Medium = Medium, nPorts = 1, use_m_flow_in = true) annotation(
          Placement(visible = true, transformation(origin = {50, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Step step(height = 0.0001, offset = 1, startTime = 1) annotation(
          Placement(visible = true, transformation(origin = {-122, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Step step1(height = -0.0001, startTime = 2) annotation(
          Placement(visible = true, transformation(origin = {-124, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Step step2(height = 0.0001, startTime = 5) annotation(
          Placement(visible = true, transformation(origin = {76, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Math.Add add annotation(
          Placement(visible = true, transformation(origin = {-86, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
          Placement(visible = true, transformation(origin = {2, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Components.Mass mass(m = 1) annotation(
          Placement(visible = true, transformation(origin = {-36, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Components.Mass mass1(m = 1) annotation(
          Placement(visible = true, transformation(origin = {40, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(boundary.ports[1], linearActuatorDynamic1.port_a) annotation(
          Line(points = {{-40, 30}, {0, 30}, {0, 0}}, color = {0, 127, 255}));
        connect(boundary1.ports[1], linearActuatorDynamic1.port_a) annotation(
          Line(points = {{40, 30}, {0, 30}, {0, 0}}, color = {0, 127, 255}));
        connect(step.y, add.u1) annotation(
          Line(points = {{-110, 70}, {-98, 70}, {-98, 56}}, color = {0, 0, 127}));
        connect(step1.y, add.u2) annotation(
          Line(points = {{-112, 30}, {-98, 30}, {-98, 44}}, color = {0, 0, 127}));
        connect(add.y, boundary.m_flow_in) annotation(
          Line(points = {{-74, 50}, {-66, 50}, {-66, 38}, {-60, 38}}, color = {0, 0, 127}));
        connect(step2.y, boundary1.m_flow_in) annotation(
          Line(points = {{88, 70}, {94, 70}, {94, 38}, {60, 38}}, color = {0, 0, 127}));
        connect(mass.flange_b, linearActuatorDynamic1.flange_a) annotation(
          Line(points = {{-26, -10}, {-10, -10}}, color = {0, 127, 0}));
        connect(linearActuatorDynamic1.flange_b, mass1.flange_a) annotation(
          Line(points = {{10, -10}, {30, -10}}, color = {0, 127, 0}));
        annotation(
          Diagram(coordinateSystem(extent = {{-140, 80}, {100, -20}})),
          experiment(StartTime = 0, StopTime = 6, Tolerance = 1e-06, Interval = 0.012));
      end DynamicTester;
    end MechAndFluidTests;

    model ValveControlTest
      Components.ValveControl.LUTControl lUTControl annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor[6](each R = 10) annotation(
        Placement(visible = true, transformation(origin = {0, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground[6] annotation(
        Placement(visible = true, transformation(origin = {-30, -20}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    equation
      connect(lUTControl.pin_p, resistor.n) annotation(
        Line(points = {{10, 0}, {20, 0}, {20, -40}, {10, -40}}, color = {0, 0, 255}, thickness = 0.5));
      connect(resistor.p, ground.p) annotation(
        Line(points = {{-10, -40}, {-20, -40}, {-20, -20}}, color = {0, 0, 255}));
      connect(lUTControl.pin_n, ground.p) annotation(
        Line(points = {{-10, 0}, {-20, 0}, {-20, -20}}, color = {0, 0, 255}));
      annotation(
        experiment(StartTime = 0, StopTime = 610, Tolerance = 1e-6, Interval = 1.22));
    end ValveControlTest;

    package FluidTests
      model MflowTest
        replaceable package Medium = Modelica.Media.Air.DryAirNasa;
        Modelica.Blocks.Logical.OnOffController onOffController(bandwidth = 50000) annotation(
          Placement(visible = true, transformation(origin = {24, 156}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = 0.001, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nPorts = 3, use_portsData = false) annotation(
          Placement(visible = true, transformation(origin = {-46, 116}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
          Placement(visible = true, transformation(origin = {-16, 156}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear(redeclare package Medium = Medium, dp_nominal = 49999.99999999999, m_flow_nominal = 0.01) annotation(
          Placement(visible = true, transformation(origin = {20, 96}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression p_lim(y = 175000) annotation(
          Placement(visible = true, transformation(origin = {-16, 186}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Logical.Switch switch1 annotation(
          Placement(visible = true, transformation(origin = {74, 156}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression realExpression(y = 1) annotation(
          Placement(visible = true, transformation(origin = {30, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.FixedBoundary ambient(redeclare package Medium = Medium, nPorts = 2) annotation(
          Placement(visible = true, transformation(origin = {56, 66}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression realExpression1(y = 0) annotation(
          Placement(visible = true, transformation(origin = {30, 188}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Sources.MassFlowSource_T pressurizedAir(redeclare package Medium = Medium, m_flow = 0.002, nPorts = 2) annotation(
          Placement(visible = true, transformation(origin = {-98, 96}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial) annotation(
          Placement(visible = true, transformation(origin = {-98, 184}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Fluid.Valves.ValveLinear valveLinear1(redeclare package Medium = Medium, dp_nominal = 49999.99999999999, m_flow_nominal = 0.01) annotation(
          Placement(visible = true, transformation(origin = {-50, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression realExpression2(y = 0) annotation(
          Placement(visible = true, transformation(origin = {-94, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(valveLinear.port_b, ambient.ports[1]) annotation(
          Line(points = {{30, 96}, {38, 96}, {38, 66}, {46, 66}}, color = {0, 127, 255}));
        connect(p_lim.y, onOffController.reference) annotation(
          Line(points = {{-5, 186}, {3, 186}, {3, 162}, {11, 162}}, color = {0, 0, 127}));
        connect(pressure.p, onOffController.u) annotation(
          Line(points = {{-5, 156}, {3, 156}, {3, 150}, {11, 150}}, color = {0, 0, 127}));
        connect(volume.ports[3], pressure.port) annotation(
          Line(points = {{-46, 106}, {-46, 96}, {-16, 96}, {-16, 146}}, color = {0, 127, 255}));
        connect(switch1.y, valveLinear.opening) annotation(
          Line(points = {{85, 156}, {103, 156}, {103, 116}, {19, 116}, {19, 104}}, color = {0, 0, 127}));
        connect(volume.ports[1], valveLinear.port_a) annotation(
          Line(points = {{-46, 106}, {-46, 96}, {10, 96}}, color = {0, 127, 255}));
        connect(onOffController.y, switch1.u2) annotation(
          Line(points = {{35, 156}, {61, 156}}, color = {255, 0, 255}));
        connect(pressurizedAir.ports[1], volume.ports[2]) annotation(
          Line(points = {{-88, 96}, {-46, 96}, {-46, 106}}, color = {0, 127, 255}));
        connect(realExpression1.y, switch1.u1) annotation(
          Line(points = {{42, 188}, {50, 188}, {50, 164}, {62, 164}}, color = {0, 0, 127}));
        connect(realExpression.y, switch1.u3) annotation(
          Line(points = {{42, 130}, {52, 130}, {52, 148}, {62, 148}}, color = {0, 0, 127}));
        connect(pressurizedAir.ports[2], valveLinear1.port_a) annotation(
          Line(points = {{-88, 96}, {-70, 96}, {-70, 60}, {-60, 60}}, color = {0, 127, 255}));
        connect(valveLinear1.port_b, ambient.ports[2]) annotation(
          Line(points = {{-40, 60}, {2, 60}, {2, 66}, {46, 66}}, color = {0, 127, 255}));
        connect(realExpression2.y, valveLinear1.opening) annotation(
          Line(points = {{-82, 78}, {-50, 78}, {-50, 68}}, color = {0, 0, 127}));
        annotation(
          Diagram(coordinateSystem(extent = {{-120, 200}, {100, 60}})),
          experiment(StartTime = 0, StopTime = 5, Tolerance = 1e-06, Interval = 0.002));
      end MflowTest;
    end FluidTests;
  end Tests;

  model test1
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
      Placement(visible = true, transformation(origin = {-32, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Inductor inductor(L = 10)  annotation(
      Placement(visible = true, transformation(origin = {-26, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Sources.StepVoltage stepVoltage(V = 12, startTime = 0.5)  annotation(
      Placement(visible = true, transformation(origin = {-70, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(inductor.n, ground.p) annotation(
      Line(points = {{-16, 56}, {0, 56}, {0, 16}, {-32, 16}}, color = {0, 0, 255}));
  connect(stepVoltage.n, ground.p) annotation(
      Line(points = {{-60, 28}, {-32, 28}, {-32, 16}}, color = {0, 0, 255}));
  connect(stepVoltage.p, inductor.p) annotation(
      Line(points = {{-80, 28}, {-80, 56}, {-36, 56}}, color = {0, 0, 255}));
  end test1;
  annotation(
    uses(Modelica(version = "4.0.0")));
end WormBot;
