function FA = calcFA(pulse,Amat)
    % in deg!
    FA = rad2deg(Amat * pulse(:));
end