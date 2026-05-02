function [p,q] = calc_power(state)
    
    theta_gfm = state(:,9);

    i_LineD = state(:,11);
    i_LineQ = state(:,12);

    v_od = state(:,1);
    v_oq = state(:,2);

    % Line current is exactly the GFM output current i_o
    i_od = i_LineD .* cos(theta_gfm) + i_LineQ .* sin(theta_gfm);
    i_oq = -i_LineD .* sin(theta_gfm) + i_LineQ .* cos(theta_gfm);

    % Power calculation
    p = 1.5 * (v_od .* i_od + v_oq .* i_oq);
    q = 1.5 * (v_oq .* i_od - v_od .* i_oq); % Note: q = vq*id - vd*iq
end