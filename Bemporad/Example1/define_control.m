function [ Kx, Kc ] = define_control(G, W, S, G_tio, W_tio, S_tio, H, F)

    T = inv(H)*G_tio'*inv(G_tio*inv(H)*G_tio')
    Kx = (T*S_tio-inv(H)*F');
    Kc = T*W_tio;
    Kx = Kx(1,:);
    Kc = Kc(1,:);
end

