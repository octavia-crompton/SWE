function [C,K,theta] = vanGenuchten(h,phi) 
    alpha   = phi(1);
    theta_S = phi(2);
    theta_R = phi(3);
    n       = phi(4);
    m       = phi(5);
    Ksat    = phi(6); 
    % Compute the volumetric moisture content
    theta = (theta_S -theta_R)./(1 + (alpha.*abs(h)).^n).^m + theta_R; 
    % Compute the effective saturation
    Se = ((theta - theta_R)./(theta_S - theta_R)); 
    % Compute the hydraulic conductivity
    K = Ksat.*Se.^(1/2).*(1 - (1 - Se.^(1/m)).^m).^2;
    % Compute the specific moisture storage
    C = -alpha.*n.*sign(h).*(1/n - 1).*(alpha.*abs(h)).^(n - 1).*(theta_R - ...
        theta_S).*((alpha.*abs(h)).^n + 1).^(1/n - 2);
%     for i=1:length(h)
%         if h(i) > 0
%             C(i) = 0;
%             theta(i) = (theta_S -theta_R)./(1 + (alpha.*0).^n).^m + theta_R;
%             K(i) = Ksat;
%         end
%     end
end


