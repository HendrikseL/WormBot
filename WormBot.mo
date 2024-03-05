package WormBot
   
package Components

    package Valves
      partial model ValveInterface
      end ValveInterface;
    end Valves;
    
    package SphericalActuators
      partial model SphereIterface
      end SphereIterface;
        end SphericalActuators;
    
    package CylindricalActuators
  partial model CylinderInterface
      end CylinderInterface;
    end CylindricalActuators;

    package ValveControl
  partial model ControlInterface
      end ControlInterface;
    end ValveControl;
  end Components;

  package Systems
  end Systems;

  package Tests
  end Tests;
end WormBot;
