classdef AircraftModel
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
       xcg
       W
       m
       S
       c
       b
       V
       h
       rho
       muc
       mub
       Kxx
       Kzz
       Kxz
       Kyy
       Cx0
       Cxu
       Cxa
       Cxq
       Cxd
       Cz0
       Czu
       Cza
       Czadot
       Czq
       Czd
       Cmu
       Cma
       Cmadot
       Cmq
       Cmd
       Cyb
       Cybdot
       Cyp 
       Cyr
       Cyda
       Cydr
       Clb
       Clp
       Clr
       Clda
       Cldr
       Cnb
       Cnbdot
       Cnp
       Cnr
       Cnda
       Cndr
       CL
       
    end
    
    methods(Static)
        function obj = AircraftModel()
            %PARAMETERS Construct an instance of this class
            %   Detailed explanation goes here
            obj.W = 54298; % N
            obj.m = 5535; % kg 
            obj.S = 30; % m^2
            obj.c = 2.057; % m
            obj.b = 15.762; % m
            obj.V = 125.675; % m/s
            obj.h = 500; % m
            obj.rho = 1000; % kg/m^3
            obj.muc = 89.7;
            obj.mub = 11.7;
            obj.Kxx = 0.0141; % value already squared
            obj.Kzz = 0.0378; % value already squared
            obj.Kxz = 0.0011; 
            obj.Kyy = 1.4996; % value already squared
            obj.Cx0 = 0;
            obj.Cxu = -0.0032;
            obj.Cxa = 0.1692;
            obj.Cxq = -0.0450;
            obj.Cxd = 0;
            obj.Cz0 = -0.2292;
            obj.Czu = -0.4592;
            obj.Cza = -5.7874;
            obj.Czadot = -0.3980;
            obj.Czq = -4.5499;
            obj.Czd = -0.5798;
            obj.Cmu = 0.0236;
            obj.Cma = -0.7486;
            obj.Cmadot = -4.2255;
            obj.Cmq = -7.4647;
            obj.Cmd = -1.444;
            obj.Cyb = -0.4046;
            obj.Cybdot = 0;
            obj.Cyp = -0.0733;
            obj.Cyr = 0.1193;
            obj.Cyda = 0;
            obj.Cydr = 0.3037;
            obj.Clb = -0.1090;
            obj.Clp = -0.5194;
            obj.Clr = 0.1039;
            obj.Clda = -0.2292;
            obj.Cldr = 0.0446;
            obj.Cnb = 0.06786;
            obj.Cnbdot = 0;
            obj.Cnp = 0.0010;
            obj.Cnr = -0.1279;
            obj.Cnda = 0.0071;
            obj.Cndr = -0.1261;
            obj.CL = obj.W/(1/2 * obj.rho * obj.V^2 * obj.S);
        end
    end
    methods
        function system = state_space(obj)
            C_1 = [(obj.Cybdot - 2 * obj.mub) * obj.b/obj.V, 0, 0, 0;
                   0, -1 * obj.b/(2*obj.V), 0, 0;
                   0, 0, -2 * obj.mub * obj.Kxx * obj.b^2/obj.V^2, 2 * obj.mub * obj.Kxz * obj.b^2/obj.V^2;
                   obj.Cnbdot * obj.b/obj.V, 0, 2 * obj.b^2 * obj.mub * obj.Kxz/obj.V^2,  ...
                   -2 * obj.b^2 * obj.mub * obj.Kzz/obj.V^2];

            C_2 = [obj.Cyb, obj.CL, obj.Cyp * obj.b/(2 * obj.V), ...
                             obj.b/(2 * obj.V) * (obj.Cyr - 4 * obj.mub);
                            0, 0,  obj.b/(2*obj.V), 0;
                            obj.Clb, 0, obj.b/(obj.V * 2) * obj.Clp, obj.b/(obj.V * 2) * obj.Clr;
                            obj.Cnb, 0, obj.b/(2 * obj.V) * obj.Cnp, ...
                             obj.b/(2 * obj.V) * obj.Cnr];

            C_3 = [obj.Cyda, obj.Cydr;
                            0, 0;
                            obj.Clda, obj.Cldr;
                            obj.Cnda, obj.Cndr;];
            A = -C_1\C_2;
            B = -C_1\C_3;
            C = eye(4);
            D = 0;
            system = ss(A, B, C, D, 'StateName', {'\beta'; '\phi'; 'p'; 'r'}, ... 
                                     'OutputName', {'\beta'; '\phi'; 'p'; 'r'},  ...
                                      'InputName', {'\delta_a'; '\delta_r'});
            
        end
        
        function system = simple_state_space(obj)
            a_2 = -2 * obj.b^2 * obj.mub * obj.Kzz/obj.V^2;
            a_1 = obj.b/(2 * obj.V) * obj.Cnr;
            a_0 = -obj.Cnb;
            A = [0, (2 * obj.V)/obj.b; -a_0/a_2 * obj.b/(2 * obj.V) , -a_1/a_2];
            B = [0, 0; obj.Cnda, obj.Cndr];
%             B = -C_1\C_3;
            C = [1, 0; 0, 1];
            D = 0;
            system = ss(A, B, C, D, 'StateName', {'\psi'; 'r'}, ...
                                    'OutputName',{'\psi'; 'r'}, ...
                                      'InputName', {'\delta_a'; '\delta_r'});
        end
        
        
        function system = augmented_state_space(obj)
            cit2a;
            C = zeros(4, 10);
            C(1,1) = 1; C(2,2) = 1; C(3,3) = 1; C(4,4) = 1;
            D = 0;
            system = ss(A, B, C, D);
        end
        
        function system = augmented_simple_state_space(obj)
            cit2a;
            A12 = [zeros(1, 6); A(4, 5:end)];
            A21 = zeros(6, 2);
            A22 = A(5:end, 5:end);
            system_red = obj.simple_state_space();
            A_new = [system_red.A, A12;
                    A21, A22];
            B12 = zeros(2, 3);
            B22 = B(5:end, 3:end);
            B21 = zeros(6, 2);
            B_new = [system_red.B, B12; B21, B22];
            C12 = zeros(2, 6);
            C_new = [system_red.C, C12];
            system = ss(A_new, B_new, C_new, 0);  
        end
        
        function updated_system = feedback(system, K)
            updated_system.A = system.A - system.B * K;
        end
    end
end
