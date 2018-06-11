% shear buckling code 
% unfinished, when it is finished, it should load 
% all panels and compare all panels between each rib.
TopX = nx(BIndex);
TopY = Topel.posY(BIndex);
BotY = Botel.posY(BIndex);
dX = diff(TopX);
TopdY = diff(TopY);
BotdY = diff(BotY);
TopStrDistance = sqrt(dX.^2 + TopdY.^2);
BotStrDistance = sqrt(dX.^2 + BotdY.^2);

v = 0.33;   % Pasion's Ratio
k = 9;

% Cycle Each Load Case for Skin Shear Buckling
Maxshearflow=zeros(1,12);
Maxshear=zeros(1,12);
   
  for i = 1:length(TopStrDistance)
    b = TopStrDistance(i);
    sigma_cr = (k*(pi^2)*E*skint^2)/(12*(1 - v^2)*b^2);

    if abs(min(min_z))< sigma_cr
        fprintf('top panel %2i: ', i)
        fprintf('b = %5.3f m ',TopStrDistance(i))
        fprintf('passed  ')
        fprintf('sigma_cr: %5.3e',sigma_cr)
        fos = sigma_cr/abs(min(min_z));
        fprintf('  fos: %4.2f\n',fos)
    else
        fprintf('not passed\n')
    end
  end

  for i = 1:length(BotStrDistance)
    b = BotStrDistance(i);
    sigma_cr = (k*(pi^2)*E*skint^2)/(12*(1 - v^2)*b^2);

    if abs(min(min_z))< sigma_cr
        fprintf('bottom panel %2i: ', i)
        fprintf('b = %5.3f m ',BotStrDistance(i))
        fprintf('passed  ')
        fprintf('sigma_cr: %5.3e',sigma_cr)
        fos = sigma_cr/abs(min(min_z));
        fprintf('  fos: %4.2f\n',fos)
    else
        fprintf('not passed\n')
    end
  end
