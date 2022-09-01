function sol = LieEuler(vector_field,exponential,action,p,h,sigma0,trajectory,t)
% Lie-Euler time integrator (Runge-Kutta-Munthe-Kaas method of order 1)
%
% :param vector_field: right hand side of the ODE [type: function handle]
% :param exponential: exponential map from the Lie algebra to the Lie group [type: function handle]
% :param action: Lie group action [type: function handle]
% :param p: solution at time t_n [type: float, 3x14 matrix]
% :param h: time step size [type: float]
% :param sigma0: initial value of the curve sigma on the Lie algebra [type: float, 30x1 vector]
% :param trajectory: desired trajectory [type: function handle]
% :param t: discrete time t_n [type: float]
%
% :returns: solution at time t_(n+1) [type: float, 3x14 matrix]

    sol = action(exponential(h*vector_field(sigma0,p,trajectory(t))),p);
end
