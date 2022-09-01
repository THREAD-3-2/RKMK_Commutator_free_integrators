function sol = CFree4(f,action,exponential,h,p,trajectory,t)
% Commutator-free time integrator of order 4
%
% :param f: map f from the phase space (on which the vector field is defined) to the Lie algebra [type: function handle]
% :param action: Lie group action [type: function handle]
% :param exponential: exponential map from the Lie algebra to the Lie group [type: function handle]
% :param h: time step size [type: float]
% :param p: solution at time t_n [type: float, 3x14 matrix]
% :param trajectory: desired trajectory [type: function handle]
% :param t: discrete time t_n [type float]
%
% :returns: solution at time t_(n+1) [type: float, 3x14 matrix]
        
    gAc = @(g,x) action(exponential(g),x); 
    
    Y1 = p;
    k1 = h*f(Y1,trajectory(t));
    
    Y2 = gAc(k1/2,p);
    k2 = h*f(Y2,trajectory(t+h/2));
    
    Y3 = gAc(k2/2,p);
    k3 = h*f(Y3,trajectory(t+h/2));
    
    Y4 = gAc(k3-k1/2,Y2);
    k4 = h*f(Y4,trajectory(t+h));
    
    yHalf = gAc(1/12*(3*k1+2*k2+2*k3-k4),p);
    
    sol = gAc(1/12*(-k1+2*k2+2*k3+3*k4),yHalf); 

end
