[x,nextRange] = create_initial_gen_Branch(Peak_torque, Min_torque);
        
nextRange = (Peak_torque - Min_torque)/2;
x(1) = Min_torque;
x(2) = Min_torque + nextRange;
x(3) = Peak_torque;
x(4) = Peak_torque;
x(5) = Min_torque + nextRange;
x(6) = Min_torque;


