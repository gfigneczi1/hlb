function [Poptim, U_left, dY_left, U_right, dY_right] = functional_driverModelLearning(GT_U, GT_Y ,ver)
    % finding the situations where the curvature is greater than a threshold
    % (2.5e-4 1/m)
    curveIndeces = abs(((GT_U(:,1)+GT_U(:,2)+GT_U(:,3))/3)) > (2.5e-4/1e-3);
   if (ver==0)
       %7x3 curvature gradient model
    p1 = inv(GT_U(:,curveIndeces)*GT_U(:,curveIndeces)')*GT_U(:,curveIndeces)*(GT_Y(2,curveIndeces))';
    p2 = inv(GT_U(:,curveIndeces)*GT_U(:,curveIndeces)')*GT_U(:,curveIndeces)*(GT_Y(3,curveIndeces))';
    p3 = inv(GT_U(:,curveIndeces)*GT_U(:,curveIndeces)')*GT_U(:,curveIndeces)*(GT_Y(4,curveIndeces))';
    Poptim = [p1'; p2'; p3'];    
   elseif(ver==1)
       %3x1 curvature model - diagonal
       c = polyfit(GT_U(1,curveIndeces),GT_Y(2,curveIndeces),1);
       p1 = c(1);
       c = polyfit(GT_U(2,curveIndeces),GT_Y(3,curveIndeces),1);
       p2 = c(1);
       c = polyfit(GT_U(3,curveIndeces),GT_Y(4,curveIndeces),1);
       p3 = c(1);
       Poptim = [[p1 0 0; 0 p2 0; 0 0 p3] zeros(3,4)];
   elseif(ver==2)
       %3x3 curvature model
       p1 = inv(GT_U(1:3,curveIndeces)*GT_U(1:3,curveIndeces)')*GT_U(1:3,curveIndeces)*(GT_Y(2,curveIndeces))';
        p2 = inv(GT_U(1:3,curveIndeces)*GT_U(1:3,curveIndeces)')*GT_U(1:3,curveIndeces)*(GT_Y(3,curveIndeces))';
        p3 = inv(GT_U(1:3,curveIndeces)*GT_U(1:3,curveIndeces)')*GT_U(1:3,curveIndeces)*(GT_Y(4,curveIndeces))';
        Poptim = [p1'; p2'; p3']; 
        Poptim = [Poptim zeros(3,4)];
   elseif(ver==3)
       %3x3 curvature model - upper triangle
       p1 = inv(GT_U(1:3,curveIndeces)*GT_U(1:3,curveIndeces)')*GT_U(1:3,curveIndeces)*(GT_Y(2,curveIndeces))';
        p2 = inv(GT_U(2:3,curveIndeces)*GT_U(2:3,curveIndeces)')*GT_U(2:3,curveIndeces)*(GT_Y(3,curveIndeces))';
        p3 = inv(GT_U(3,curveIndeces)*GT_U(3,curveIndeces)')*GT_U(3,curveIndeces)*(GT_Y(4,curveIndeces))';
        Poptim = [p1'; [0 p2']; [0 0 p3]]; 
        Poptim = [Poptim zeros(3,4)];
   elseif(ver==4)
       % diagonal with one curvature gradient
       p1 = inv([GT_U(1,curveIndeces); GT_U(6,curveIndeces)]*[GT_U(1,curveIndeces); GT_U(6,curveIndeces)]')*[GT_U(1,curveIndeces); GT_U(6,curveIndeces)]*(GT_Y(2,curveIndeces))';
       p2 = inv(GT_U(2,curveIndeces)*GT_U(2,curveIndeces)')*GT_U(2,curveIndeces)*(GT_Y(3,curveIndeces))';
       p3 = inv(GT_U(3,curveIndeces)*GT_U(3,curveIndeces)')*GT_U(3,curveIndeces)*(GT_Y(4,curveIndeces))';
       Poptim = [[p1(1) 0 0 0 0 p1(2) 0]; [0 p2 0 0 0 0 0]; [0 0 p3 0 0 0 0]];
   elseif (ver==5)
        % 3x3 model with one curvature gradient
        p1 = inv([GT_U(1:3,curveIndeces); GT_U(6,curveIndeces)]*[GT_U(1:3,curveIndeces); GT_U(6,curveIndeces)]')*[GT_U(1:3,curveIndeces); GT_U(6,curveIndeces)]*(GT_Y(2,curveIndeces))';
        p2 = inv(GT_U(1:3,curveIndeces)*GT_U(1:3,curveIndeces)')*GT_U(1:3,curveIndeces)*(GT_Y(3,curveIndeces))';
        p3 = inv(GT_U(1:3,curveIndeces)*GT_U(1:3,curveIndeces)')*GT_U(1:3,curveIndeces)*(GT_Y(4,curveIndeces))';
        Poptim = [[p1(1:3)' 0 0 p1(4) 0]; [p2' 0 0 0 0]; [p3' 0 0 0 0]];
   elseif (ver==6)
        % upper triangle model with one curvture gradient
        p1 = inv([GT_U(1:3,curveIndeces); GT_U(6,curveIndeces)]*[GT_U(1:3,curveIndeces); GT_U(6,curveIndeces)]')*[GT_U(1:3,curveIndeces); GT_U(6,curveIndeces)]*(GT_Y(2,curveIndeces))';
        p2 = inv(GT_U(2:3,curveIndeces)*GT_U(2:3,curveIndeces)')*GT_U(2:3,curveIndeces)*(GT_Y(3,curveIndeces))';
        p3 = inv(GT_U(3,curveIndeces)*GT_U(3,curveIndeces)')*GT_U(3,curveIndeces)*(GT_Y(4,curveIndeces))';
        Poptim = [[p1(1:3)' 0 0 p1(4) 0]; [0 p2' 0 0 0 0]; [0 0 p3 0 0 0 0]];
   elseif (ver==7)
       % E-LDM
       positiveIndeces = (GT_U(:,1)+GT_U(:,2)+GT_U(:,3))/3 > (2.5e-4/1e-3);
       negativeIndeces = (GT_U(:,1)+GT_U(:,2)+GT_U(:,3))/3 < -(2.5e-4/1e-3);

       if (isempty(find(curveIndeces==0)))
            p0 = 0;
       else
            p0 = mean(mean(GT_Y(curveIndeces==0,:)));
       end
       dY = GT_Y-p0;
       if (isempty(find(positiveIndeces>0)))
           pLeft = zeros(3,3);
           U_left = []; dY_left = [];
       else
           U_left = GT_U(positiveIndeces,:);
           dY_left = GT_Y(positiveIndeces,:);
           p1 = inv(U_left'*U_left)*U_left'*dY_left(:,1);
           p2 = inv(U_left'*U_left)*U_left'*dY_left(:,2);
           p3 = inv(U_left'*U_left)*U_left'*dY_left(:,3);
           pLeft = [p1 p2 p3];
       end
       if (isempty(find(negativeIndeces>0)))
           pRight = zeros(3,3);
           U_right = []; dY_right = [];
       else
           U_right = GT_U(negativeIndeces,:);
           dY_right = GT_Y(negativeIndeces,:);
           p4 = inv(U_right'*U_right)*U_right'*dY_right(:,1);
           p5 = inv(U_right'*U_right)*U_right'*dY_right(:,2);
           p6 = inv(U_right'*U_right)*U_right'*dY_right(:,3);
           pRight = [p4 p5 p6];
       end
       Poptim = [reshape(pLeft,9,1); reshape(pRight,9,1); p0];
   elseif (ver==8)
       % E-LDM from pure inputs
       positiveIndeces = (GT_U(:,1)+GT_U(:,2)+GT_U(:,3))/3 > (2.5e-4/1e-3);
       negativeIndeces = (GT_U(:,1)+GT_U(:,2)+GT_U(:,3))/3 < -(2.5e-4/1e-3);

       if (isempty(find(curveIndeces==0)))
            p0 = 0;
       else
            p0 = mean(mean(GT_Y(curveIndeces==0,:)));
       end
       dY = GT_Y-p0;
       if (isempty(find(positiveIndeces>0)))
           pLeft = zeros(3,3);
           U_left = []; dY_left = [];
       else
           U_left = GT_U(positiveIndeces,:);
           dY_left = GT_Y(positiveIndeces,:);
           pLeft = inv(U_left'*U_left)*U_left'*dY_left;
       end
       if (isempty(find(negativeIndeces>0)))
           pRight = zeros(3,3);
           U_right = []; dY_right = [];
       else
           U_right = GT_U(negativeIndeces,:);
           dY_right = GT_Y(negativeIndeces,:);
           pRight = inv(U_right'*U_right)*U_right'*dY_right;
       end
%        p1 = inv(GT_U(1:3,positiveIndeces)*GT_U(1:3,positiveIndeces)')*GT_U(1:3,positiveIndeces)*(GT_Y(1,positiveIndeces)-p0)';
%        p2 = inv(GT_U(1:3,positiveIndeces)*GT_U(1:3,positiveIndeces)')*GT_U(1:3,positiveIndeces)*(GT_Y(2,positiveIndeces)-p0)';
%        p3 = inv(GT_U(1:3,positiveIndeces)*GT_U(1:3,positiveIndeces)')*GT_U(1:3,positiveIndeces)*(GT_Y(3,positiveIndeces)-p0)';
%        
%        p4 = inv(GT_U(1:3,negativeIndeces)*GT_U(1:3,negativeIndeces)')*GT_U(1:3,negativeIndeces)*(GT_Y(1,negativeIndeces)-p0)';
%        p5 = inv(GT_U(1:3,negativeIndeces)*GT_U(1:3,negativeIndeces)')*GT_U(1:3,negativeIndeces)*(GT_Y(2,negativeIndeces)-p0)';
%        p6 = inv(GT_U(1:3,negativeIndeces)*GT_U(1:3,negativeIndeces)')*GT_U(1:3,negativeIndeces)*(GT_Y(3,negativeIndeces)-p0)';
       
       Poptim = [reshape(pLeft,9,1); reshape(pRight,9,1); p0]; % p4' p0; p2' p5' p0; p3' p6' p0];
   end
   Poptim = reshape(Poptim,size(Poptim,1)*size(Poptim,2),1);
end

