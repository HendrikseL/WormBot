package WormBot
  package Components
    package Valves
      partial model ValveInterface
      end ValveInterface;
    end Valves;

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
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium model";
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
    end CylindricalActuators;

    package ValveControl
      partial model ControlInterface
      end ControlInterface;
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
        Modelica.Blocks.Sources.RealExpression nearZero(y = 10^(-6)) annotation(
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
    end Friction;
  end Components;

  package Systems
  end Systems;

  package Tests
    package MechanicalTests
      model FrictionAndOnOffTest
        parameter Modelica.Units.SI.Mass m = 1;
        Modelica.Mechanics.Translational.Sources.Force force annotation(
          Placement(visible = true, transformation(origin = {60, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Ramp ramp(duration = 0, height = 3, offset = 0, startTime = 0) annotation(
          Placement(visible = true, transformation(origin = {70, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.LuGre frictionActuator annotation(
          Placement(visible = true, transformation(origin = {20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Components.Mass mass(m = m) annotation(
          Placement(visible = true, transformation(origin = {-50, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression realExpression(y = m*9.81) annotation(
          Placement(visible = true, transformation(origin = {-80, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        WormBot.Components.Friction.LuGre frictionJoint(mu_k = 0.2, mu_s = 0.25) annotation(
          Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.BooleanStep booleanStep(startTime = 1) annotation(
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
          experiment(StartTime = 0, StopTime = 3, Tolerance = 1e-06, Interval = 1e-05));
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
    end MechAndFluidTests;
  end Tests;
  annotation(
    uses(Modelica(version = "4.0.0")));
end WormBot;
