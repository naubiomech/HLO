function varargout = Optimization(varargin)
% OPTIMIZATION MATLAB code for Optimization.fig
%      OPTIMIZATION, by itself, creates a new OPTIMIZATION or raises the existing
%      singleton*.
%
%      H = OPTIMIZATION returns the handle to a new OPTIMIZATION or the handle to
%      the existing singleton*.
%
%      OPTIMIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPTIMIZATION.M with the given input arguments.
%
%      OPTIMIZATION('Property','Value',...) creates a new OPTIMIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Optimization_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Optimization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Optimization

% Last Modified by GUIDE v2.5 27-Mar-2019 15:08:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Optimization_OpeningFcn, ...
                   'gui_OutputFcn',  @Optimization_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end


% --- Executes just before Optimization is made visible.
function Optimization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Optimization (see VARARGIN)

% Choose default command line output for Optimization
handles.output = hObject;

close(figure(100),figure(101),figure(102));

global GUI_Variables

disp('Initializing optimization GUI...')

tTCP =  tcpip('10.18.48.77',55555,'networkrole','client',...
    'InputBufferSize',4096,'Timeout',60); % Omnia TCP
eTCP = tcpip('0.0.0.0',3000,'networkrole','server');       % A_EXO TCP
enablestream = '<OmniaXB><System><EnableRealTimeInformation><Enabled>1</Enabled></EnableRealTimeInformation></System></OmniaXB>';

GUI_Variables = struct('OmniaTCP',tTCP,'EXOTCP',eTCP,'StreamStr',enablestream,...
    'GenNum',1,'MidGen',0,'CompleteCond',0,'SubjectMass',1,'PkTRQ',0.35,'MinTRQ',0.2,...
    'SSID',NaN,'NumParams',0,'CondPerGen',7,'ConditionTime',240,'Stopped',0,...
    'TestDate',' ','pullDir',' ','SeedNextGen',0,'Streaming',0);

set(handles.GenNumber,'String','1');
set(handles.MidGenCheckbox,'Value',0);
set(handles.LastConditionCompleted,'enable','off');
set(handles.LastConditionCompleted,'String','0');
set(handles.ConditionTime,'String','240');
set(handles.SubjectMass,'String','1');
set(handles.MaxPeakTorque,'String','0.35');
set(handles.MinTRQ,'String','0.2');
set(handles.SubjectID,'String','SSI_');
set(handles.OmniaTCPStatus,'Color',[1 1 1]);
set(handles.ExoTCPStatus,'Color',[1 1 1]);
set(handles.OneParamCheck,'Value',0);
set(handles.TwoParamCheck,'Value',0);
set(handles.TestDate,'enable','off');
set(handles.TestDate,'String','DD-Mon-YYYY');
set(handles.StatusText,'String',...
    ['Initialization complete. Ready to begin. '...
    'Input optimizer and subject settings, open TCP connections, '...
    'and click "Start Optimizer" button.'])

disp('Initialization complete.')
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Optimization wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = Optimization_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end


% --- Executes on button press in MidGenCheckbox.
function MidGenCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to MidGenCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

try
    GUI_Variables.MidGen = get(hObject,'Value');
    if GUI_Variables.MidGen == 1
        set(handles.LastConditionCompleted,'enable','on');
    else
        set(handles.LastConditionCompleted,'enable','off');
        set(handles.LastConditionCompleted,'string','0');
        GUI_Variables.CompleteCond = 0;
    end
catch
end

if ~any(isnan(GUI_Variables.SSID)) && (GUI_Variables.GenNum>1 || GUI_Variables.MidGen==1)
    set(handles.TestDate,'enable','on')
elseif GUI_Variables.GenNum == 1 && GUI_Variables.MidGen == 0
    set(handles.TestDate,'enable','off')
    set(handles.TestDate,'String','DD-Mon-YYYY')
    GUI_Variables.pullDir = ' ';
    GUI_Variables.TestDate = ' ';
end

end

% Hint: get(hObject,'Value') returns toggle state of MidGenCheckbox


% --- Executes on button press in ExoOpenTCP.
function ExoOpenTCP_Callback(hObject, eventdata, handles)
% hObject    handle to ExoOpenTCP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

try
    eTCP = GUI_Variables.EXOTCP;
    disp('Attempting to make a TCP/IP connection to the client')
    set(handles.StatusText,'String','Attempting to make TCP connection to A_EXO...')
    pause(.01);
    fopen(eTCP);
catch
    set(handles.ExoTCPStatus,'Color',[1 0 0]);
    set(handles.StatusText,'String','TCP connection to A_EXO failed. Check client.');
end

if exist('eTCP','var')
    if eTCP.Status == "open"
        set(handles.ExoTCPStatus,'Color',[0 1 0]);
        disp('Successfully connected to A_EXO computer')
        set(handles.StatusText,'String','Successfully connected to A_EXO')
    elseif eTCP.Status == "closed"
        set(handles.ExoTCPStatus,'Color',[1 0 0]);
        set(handles.StatusText,'String','TCP connection to A_EXO failed. Check client.');
    end
end
end

% --- Executes on button press in ExoCloseTCP.
function ExoCloseTCP_Callback(hObject, eventdata, handles)
% hObject    handle to ExoCloseTCP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

try
    eTCP = GUI_Variables.EXOTCP;
catch
    set(handles.StatusText,'String','No TCP connection for A_EXO defined.');
end

if exist('eTCP','var')
    if eTCP.Status == "open"
        fclose(eTCP);
        set(handles.ExoTCPStatus,'Color',[1 1 1]);
        disp('Closed A_EXO TCP connection.');
        set(handles.StatusText,'String','Closed A_EXO TCP connection.');
    elseif eTCP.Status == "closed"
        set(handles.StatusText,'String','No A_EXO TCP connection to close.');
    end
end
end


% --- Executes on button press in OmniaOpenTCP.
function OmniaOpenTCP_Callback(hObject, eventdata, handles)
% hObject    handle to OmniaOpenTCP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

try
    t = GUI_Variables.OmniaTCP;
    str = GUI_Variables.StreamStr;
    disp('Attempting to make a TCP/IP connection to Omnia host')
    set(handles.StatusText,'String',...
        'Attempting to make TCP connection to Omnia...')
    pause(.01);
    fopen(t);
    fwrite(t,str);
catch
    set(handles.OmniaTCPStatus,'Color',[1 0 0]);
    set(handles.StatusText,'String',...
        'TCP connection to Omnia failed. Check that Omnia is open.');
end

if exist('t','var')
    if t.Status == "open"
        set(handles.OmniaTCPStatus,'Color',[0 1 0]);
        disp('Successfully connected to Omnia')
        set(handles.StatusText,'String','Successfully connected to Omnia')
    elseif t.Status == "closed"
        set(handles.OmniaTCPStatus,'Color',[1 0 0]);
        set(handles.StatusText,'String',...
            'TCP connection to Omnia failed. Check that Omnia is open.');
    end
end
end

% --- Executes on button press in OmniaCloseTCP.
function OmniaCloseTCP_Callback(hObject, eventdata, handles)
% hObject    handle to OmniaCloseTCP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

try
    t = GUI_Variables.OmniaTCP;
catch
    set(handles.StatusText,'String','No TCP connection for Omnia defined.');
end

if exist('t','var')
    if t.Status == "open"
        fclose(t);
        set(handles.OmniaTCPStatus,'Color',[1 1 1]);
        disp('Closed Omnia TCP connection.');
        set(handles.StatusText,'String','Closed Omnia TCP connection.');
    elseif t.Status == "closed"
        set(handles.StatusText,'String','No Omnia TCP connection to close.');
    end
end
end

% --- Executes on button press in StartOptimizer.
function StartOptimizer_Callback(hObject, eventdata, handles)
% hObject    handle to StartOptimizer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

GUI_Variables.Stopped=0;

set(handles.StatusText,'String','Initializing optimization...');
pause(1);

set(handles.TRQ_SetpointText,'String',' ');
set(handles.RT_SetpointText,'String',' ');
set(handles.TRQSetpointHistoryText,'String',' ');                             
set(handles.RTSetpointHistoryText,'String',' ');

if isnan(GUI_Variables.SSID)
    set(handles.StatusText,'String',...
        'Invalid subject ID defined. Cannot begin optimization.')
% elseif GUI_Variables.OmniaTCP.Status == "closed"
%     set(handles.StatusText,'String',...
%         'Omnia TCP connection not established. Cannot begin optimization.')
elseif GUI_Variables.EXOTCP.Status == "closed"
    set(handles.StatusText,'String',...
        'A_EXO TCP connection not established. Cannot begin optimization.')
elseif GUI_Variables.NumParams == 0
    set(handles.StatusText,'String',...
        'Number of parameters to optimize not selected. Cannot begin optimization.');
else

    GUI_Variables.SeedNextGen = 0;
    GenerationNumber = GUI_Variables.GenNum;
    StartingFromMidGen = GUI_Variables.MidGen;
    NextConditionMidGen = GUI_Variables.CompleteCond;
    Subject_mass = GUI_Variables.SubjectMass;
    Peak_torque = GUI_Variables.PkTRQ*Subject_mass;
    Min_torque = GUI_Variables.MinTRQ*Subject_mass;
    SSID = GUI_Variables.SSID;
    NumberofParams = GUI_Variables.NumParams;
    ConditionsPerGen = GUI_Variables.CondPerGen;
    ConditionTime = GUI_Variables.ConditionTime;
    eTCP = GUI_Variables.EXOTCP
    tTCP =  GUI_Variables.OmniaTCP
    currentDir = cd;                                % Current directory
    saveDir = [currentDir,'\',SSID,'_',date,'\'];   % Save directory
    
    if NumberofParams == 1
        InitialGuessofParams = [0.25*Subject_mass]; %[randi([0,14]); randi([30,70]) ]; %EDIT %Must be a column of parameters %For testing %Best guess (0 to Peak torque) (0 to 100)
        InitialGuessofParams(1) = InitialGuessofParams(1)/Peak_torque; %This makes it range from 0 to 1.
    elseif NumberofParams == 2
        InitialGuessofParams = [0.25*Subject_mass; 50]; %[randi([0,14]); randi([30,70]) ]; %EDIT %Must be a column of parameters %For testing %Best guess (0 to Peak torque) (0 to 100)
        InitialGuessofParams(1) = InitialGuessofParams(1)/Peak_torque; %This makes it range from 0 to 1.
        InitialGuessofParams(2) = InitialGuessofParams(2)/100; %This makes the percentage range from 0 to 1.
    end
    sigma_start = 0.3; %Coordinate wise standard deviation (step size) %For testing. 
    pc_init = zeros(NumberofParams,1); %For testing
    ps_init = zeros(NumberofParams,1); %For testing. 
    B_init = eye(NumberofParams,NumberofParams);                 
    D_init = ones(NumberofParams,1);                     
    C_init = B_init * diag(D_init.^2) * B_init';    
    counteval = 0; %This is how many total conditions have been tested.
    stopeval = 10*ConditionsPerGen; %Total number of conditions to test, currently assumes 10 generations. 
    
    mkdir(saveDir);     % Generate a save directory on first generation

    %% If just starting so need to create first generation
    if GenerationNumber == 1 && StartingFromMidGen == 0      
        
        set(handles.StatusText,'String',...
            'Starting from first condition of generation 1');
        
        xmean_history = [];
        if NumberofParams == 1
            set(handles.XMeanTrq,'String',' ');
        elseif NumberofParams == 2
            set(handles.XMeanTrq,'String',' ');
            set(handles.XMeanRT,'String',' ');
        end
        
        if NumberofParams == 1
            Initial_Gen = @create_initial_generation_newCMA_1Param;
        elseif NumberofParams == 2
            Initial_Gen = @create_initial_generation_newCMA_2Param;
        end
        
       [N, xmean, sigma, lambda, pc, ps, B, D, C, mu, weights, mueff, ...
        invsqrtC, eigeneval, chiN, cc, cs, c1, cmu, damps, x] = ...
        Initial_Gen(...
            NumberofParams, InitialGuessofParams, sigma_start, ConditionsPerGen,...
            pc_init, ps_init, B_init, D_init, C_init, Peak_torque, Min_torque);

       %X comes out as 
      % x=[p1a p1b p1c p1d
      %     p2a p2b p2c p2d
      %     p3a p3b p3c p3d] where letter is controller and number is which parameter.
      %   So transpose x to create Params full in the format I wrote the code for. 
      Params_full=x';

      %When params_full first comes out of the CMA, it comes out as [0 1] [0 1]

      %Turn them back into actual values with peak torques and percents.
      
      if NumberofParams == 1
        Params_full(:,1) = Params_full(:,1)*Peak_torque; %Peak Torque (Nm)
      elseif NumberofParams == 2
          Params_full(:,1) = Params_full(:,1)*Peak_torque; %Peak Torque (Nm)
          Params_full(:,2) = Params_full(:,2)*100;  %Rise time (% of stance time)
      end
      
      disp('Your controllers for Gen 1:')
      disp('[Peak Torque, Rise time (%)]')  
      set(handles.StatusText,'String',...
        {'Your controllers for Gen 1: '...
         '[Peak Torque, Rise Time (%)]'...
         num2str(Params_full)});
      Params_full
       % Params_full=[P1a P2a P3a
       %             P1b P2b P3b
       %            P1c P2c P3c
       %           ... ... ...
       %          P1h P2h P3h];                 
    end

    %% If you are picking up from where you left off but did complete a full generation
    if GenerationNumber>1 && StartingFromMidGen == 0
        set(handles.StatusText,'String',...
            ['Starting from first condition of generation ',num2str(GenerationNumber)]);
        GenerationAcc = GenerationNumber;
        
        currentDir = cd;                                % Current directory
        saveDir = [currentDir,'\',SSID,'_',date,'\'];   % Save directory
        
        if strcmp(GUI_Variables.TestDate,date) || strcmp(GUI_Variables.TestDate,' ')  %Compare current date to entered date
            load(fullfile(saveDir,['Completion_of_Gen_' num2str(GenerationNumber-1),'_',SSID]));
        else
            load(fullfile(GUI_Variables.pullDir,['Completion_of_Gen_' num2str(GenerationNumber-1),'_',SSID]));
        end
        
        currentDir = cd;                                % Current directory
        saveDir = [currentDir,'\',SSID,'_',date,'\'];   % Save directory
            
        GenerationNumber = GenerationAcc;
        %This loads in all the variables from the last completed generation. 
        disp(['Your controllers for Gen ' num2str(GenerationNumber)])
        disp('[Peak Torque, Rise time (%)]')
        Params_full
        set(handles.StatusText,'String',...
            {['Your controllers for Gen ' num2str(GenerationNumber)]...
             '[Peak Torque, Rise Time (%)]'...
             num2str(Params_full)});
         
        set(handles.TRQ_SetpointText,'String',' ');
        set(handles.RT_SetpointText,'String',' ');
        set(handles.TRQSetpointHistoryText,'String',' ');
        set(handles.RTSetpointHistoryText,'String', ' ');
        if NumberofParams == 1
        set(handles.XMeanTrq,'String',...
            num2str(xmean_history(:,1)));
        elseif NumberofParams == 2
            set(handles.XMeanTrq,'String',...
                num2str(xmean_history(:,1)));
            set(handles.XMeanRT,'String',...
               num2str(xmean_history(:,2)));
        end
        
        pause(1);

    end

    %Initializations if not starting partway through a generation (applies to
    %both starting at the very beggining and when starting a new generation)
    ConditionNumber = 1;
    fullrate = zeros(1,200);
    breathtimes = zeros(1,200);
    SSdata = zeros(ConditionsPerGen, (2+NumberofParams));
    Full_Metabolic_Data_to_Save = {};

    %% If you are starting from partway through a generation
    if StartingFromMidGen==1
        
    currentDir = cd;                                % Current directory
    saveDir = [currentDir,'\',SSID,'_',date,'\'];   % Save directory    

    set(handles.StatusText,'String',...
        ['Starting from mid generation on condition ' num2str(NextConditionMidGen) ...
        ' of generation ' num2str(GenerationNumber)]);

    if strcmp(GUI_Variables.TestDate,date) || strcmp(GUI_Variables.TestDate,' ')  %Compare current date to entered date
        load(fullfile(saveDir,['All_Saved_Data_Following_Gen_', num2str(GenerationNumber), '_Cond_', num2str(NextConditionMidGen),'_',SSID]))
    else
        load(fullfile(GUI_Variables.pullDir,['All_Saved_Data_Following_Gen_', num2str(GenerationNumber), '_Cond_', num2str(NextConditionMidGen),'_',SSID]))
    end
    
    currentDir = cd;                                % Current directory
    saveDir = [currentDir,'\',SSID,'_',date,'\'];   % Save directory
    
    set(handles.TRQ_SetpointText,'String',' ');
    set(handles.RT_SetpointText,'String',' ');
    set(handles.TRQSetpointHistoryText,'String',' ');
    set(handles.RTSetpointHistoryText,'String', ' ');
    if NumberofParams == 1
        if isempty(xmean_history)
            set(handles.XMeanTrq,'String',' ');
        else
            set(handles.XMeanTrq,'String',...
                num2str(xmean_history(:,1)));
        end
    elseif NumberofParams == 2
        if isempty(xmean_history)
            set(handles.XMeanTrq,'String',' ');
            set(handles.XMeanRT,'String',' ');
        else
            set(handles.XMeanTrq,'String',...
                num2str(xmean_history(:,1)));
            set(handles.XMeanRT,'String',...
               num2str(xmean_history(:,2)));
        end
    end     
    
    end

    ParamsForCondition = Params_full(ConditionNumber, :); %First row of params_full
    ParamsRemaining = Params_full(ConditionNumber:ConditionsPerGen, :)

    set(handles.StatusText,'String',...
        {['Your remaining controllers for Gen ' num2str(GenerationNumber)]...
             '[Peak Torque, Rise Time (%)]'...
             num2str(ParamsRemaining)});

    pause(1);   
    
    %% Proceed with the optimization

    disp('Setup complete, press "Play" in COSMED and press "Enter" in MATLAB')
    set(handles.StatusText,'String',...
        'Setup complete. Press "Play" in COSMED and press "ENTER" in MATLAB');
    pause

%     working=0; %Check if working
%     %Sometimes it just doesn't start right away and you just need to try again
%     %to establish the tcpclient. So this does that until it is working.
%     set(handles.StatusText,'String',...
%         'Checking for breath data from Omnia...');
%     flushinput(tTCP);
%     while working==0
%         try
%             fopen(tTCP);
%         catch
%         end
%         pause(2)
%         if tTCP.BytesAvailable>0
%             t1=fread(tTCP,tTCP.BytesAvailable);
%             if length(t1)>0
%                 working=1;
%             end
%         else
%         end
%     end

    disp('Ready to begin optimization. Have participant begin walking and press "Enter" in MATLAB');
    set(handles.StatusText,'String',...
        'Ready to begin optimization. Have participant begin walking and press "Enter" in MATLAB');
    pause;

    %Creating the timer. 
    TimerVar = timer;
    set(TimerVar, 'executionMode', 'fixedRate')
    set(TimerVar, 'TimerFcn', {@updater})
    if NumberofParams == 1
        set(TimerVar, 'UserData', [ConditionNumber-1; 0; Params_full])
    elseif NumberofParams == 2
        set(TimerVar, 'UserData', [ConditionNumber-1 ConditionNumber-1; 0 0; Params_full])
    end
    set(TimerVar, 'StartDelay', 10) %Start delay is 10 seconds, just so you know your exo code will catch the first change.
    %UserData is in the following format: 
    %     [i i; controller you are checking; controller1; controller2; etc]
    %The external code needs to check the value of 'TimerVar.UserData(2, :)
    set(TimerVar, 'Period', ConditionTime) %Period is 120 seconds, ie 2 minutes of walking
    
    if NumberofParams == 1
        NextParams = [0];
    elseif NumberofParams == 2
        NextParams=[0 0]; %Initialization of the Nextparams variable, to be checked by exo controller code. 
    end
    
     set(handles.StatusText,'String',...
         'Waiting for start command from A_EXO...');
    pause(0.1);
    SendValueFlag = 0;
    StopFlag = ' ';
    ActiveFlag = 0;
%    start(TimerVar)
    while 1 == 1
        if eTCP.Status == "open"
            [SendValueFlag,StopFlag,ActiveFlag] = CheckValue(eTCP,StopFlag,ActiveFlag);
            if ActiveFlag == 1 %Ensures exo has started checking for optimization values
                set(handles.StatusText,'String','Received start command from A_EXO.');
                GUI_Variables.Streaming = 1;
                pause(1)
                start(TimerVar)
                break
            else
            end
        else
        end
    end

    i = 0; %Initialize the counter for storing breaths
    conditiondone = 0; %Initialize to make sure it knows condition isn't done yet.
    generationdone = 0; %Initialize to say that the generation isn't done yet. 
   if tTCP.BytesAvailable>0
       t3 = fread(tTCP,tTCP.BytesAvailable); %Just to clear it.
   end
    tic %supplementary timer for tracking breath times
    CurrentCondition = ConditionNumber;
    TimerConditionBig = TimerVar.UserData;
    TimerCondition = TimerConditionBig(1,1);

    while TimerCondition<ConditionNumber %(so until the first 10 second delay has passed)
        TimerConditionBig = TimerVar.UserData;
        TimerCondition = TimerConditionBig(1,1); %Just keep checking to see if startDelay is done
    end

    set(handles.StatusText,'String',...
        'Optimizing...');

    DoneFirstValue = 0;
    SetpointHistory = [];
    while 1 == 1
        if 1 == 1
             [SendValueFlag,StopFlag,ActiveFlag] = CheckValue(eTCP,StopFlag,ActiveFlag);
            if ActiveFlag == 1 && generationdone == 0
                while 1 == 1
                     [SendValueFlag,StopFlag,ActiveFlag] = CheckValue(eTCP,StopFlag,ActiveFlag);
                    if StopFlag == 0 && DoneFirstValue == 0 && GUI_Variables.Stopped == 0
                        NextParamsBig = TimerVar.UserData;
                        NextParams = NextParamsBig(2,:)
                        if NumberofParams == 1
                             fwrite(eTCP,num2str(NextParams(1)));
                        elseif NumberofParams == 2
                             fwrite(eTCP,horzcat(num2str(NextParams(1)),'_',num2str(NextParams(2))));
                        end
                        DoneFirstValue = 1;
                        set(handles.StatusText,'String',...
                            {'Sent setpoint to A_EXO: ' ...
                            '[Peak Torque, Rise Time (%)]'...
                             num2str(NextParams)});
                         SetpointHistory = [SetpointHistory; NextParams];
                         if NumberofParams == 1
                             set(handles.TRQ_SetpointText,'String',...
                             num2str(NextParams(1)));
                             set(handles.TRQSetpointHistoryText,'String',...
                                num2str(SetpointHistory(:,1)));
                         elseif NumberofParams == 2
                             set(handles.RT_SetpointText,'String',...
                                num2str(NextParams(2)));
                             set(handles.TRQSetpointHistoryText,'String',...
                                num2str(SetpointHistory(:,1)));
                             set(handles.RTSetpointHistoryText,'String',...
                                num2str(SetpointHistory(:,2)));
                         end
                    elseif StopFlag == 0 && DoneFirstValue == 1
                        %% Big while loop for if generation is done
                        while generationdone == 0

                            %% Smaller while loop for if condition is done.
                            while conditiondone == 0
                                %every .1 seconds, run the run_on_timer function. 
                                CurrentConditionBig = TimerVar.UserData;
                                CurrentCondition = CurrentConditionBig(1,1);
                                if ConditionNumber<CurrentCondition
                                    conditiondone = 1; %Once the timer has fired so it updated the ConditionNumber.
                                else %Else you collect breath data.
                                    [tTCP,i,fullrate,breathtimes] = run_on_timer(tTCP,i, fullrate, breathtimes);
                                    pause(.1) %You will stay in this loop until you are done with the condition. 
                                    %This just searches for new breath data and stores it if there is
                                    %any. 
                                    if toc>=ConditionTime && CurrentCondition == 1 % 10 second timer offset for first condition
                                        set(handles.StatusText,'String',...
                                            'New condition in 10 seconds!');
                                    elseif toc>=ConditionTime-5 && CurrentCondition > 1
                                        set(handles.StatusText,'String',...
                                            'New condition in 5 seconds!');
                                    end
                                end
                            end
                            %Once the condition is done:
                            NextParamsBig = TimerVar.UserData;
                            NextParams = NextParamsBig(2,:)
                            if NumberofParams == 1
                                 fwrite(eTCP,num2str(NextParams(1)))
                            elseif NumberofParams == 2
                                 fwrite(eTCP,horzcat(num2str(NextParams(1)),'_',num2str(NextParams(2))));
                            end
                            set(handles.StatusText,'String',...
                            {'Sent setpoint to A_EXO: ' ...
                            '[Peak Torque, Rise Time (%)]'...
                             num2str(NextParams),...
                             ['Controller changing to next setpoint. '....
                             'To stop mid-generation after completion of '...
                             'current condition press "Pause Optimizer"']...
                             ['To seed next generation complete three conditions '...
                             'and press "Seed Next Gen"']});
                             if NumberofParams == 1
                                set(handles.TRQ_SetpointText,'String',...
                                    num2str(NextParams(1)));
                             elseif NumberofParams == 2
                                 set(handles.TRQ_SetpointText,'String',...
                                    num2str(NextParams(1)));
                                 set(handles.RT_SetpointText,'String',...
                                    num2str(NextParams(2)));
                             end
                             SetpointHistory = [SetpointHistory; NextParams];
                             if NumberofParams == 1
                                 set(handles.TRQSetpointHistoryText,'String',...
                                 num2str(SetpointHistory(:,1)));
                             elseif NumberofParams == 2
                                 set(handles.TRQSetpointHistoryText,'String',...
                                 num2str(SetpointHistory(:,1)));
                                 set(handles.RTSetpointHistoryText,'String',...
                                 num2str(SetpointHistory(:,2)));
                             end
                            tic %reinitialize timer once you start next controller
                            %nonzeros(breathtimes) %uncomment if you want to see that its working
                            [SS, y_bar]=fake_metabolic_fit_JB(ParamsForCondition); %Used for
                            Full_y_bar{ConditionNumber} = y_bar';
                            %testing the optimization
                            Full_Metabolic_Data_to_Save{ConditionNumber} = [y_bar',10*y_bar'];
                            SSdata(ConditionNumber, :) = [SS, ConditionNumber, ParamsForCondition]; 
                            ConditionNumber = ConditionNumber+1;
                            counteval = counteval+1; %Keeps track of how many conditions have been tested. 
                            %Automatically save it after every condition, so you don't have to redo
                            %it if you stop mid generation. You can delete the extra saved
                            %condition date later.
                            set(handles.StatusText,'String',...
                            {'Sent setpoint to A_EXO: ' ...
                            '[Peak Torque, Rise Time (%)]'...
                             num2str(NextParams),...
                             ['Controller changing to next setpoint. '....
                             'To stop mid-generation after completion of '...
                             'current condition press "Pause Optimizer"']...
                             ['To seed next generation complete three conditions '...
                             'and press "Seed Next Gen"'],...
                             ['Saving file "All_Saved_Data_Following_Gen_',...
                             num2str(GenerationNumber), '_Cond_',...
                                num2str(ConditionNumber-1),'_',SSID,'"']});
                            VarsToIgnore = ['eventdata|GUI_Variables|handles|hObject|ActiveFlag'...
                            '|DoneFirstValue|SendValueFlag|SetpointHistory|eTCP|tTCP'];
                            save(fullfile(saveDir,['All_Saved_Data_Following_Gen_', num2str(GenerationNumber),...
                                '_Cond_', num2str(ConditionNumber-1),'_',SSID]),...
                                '-regexp',['^(?!',VarsToIgnore,'$).']);
                             [SendValueFlag,StopFlag,ActiveFlag] = CheckValue(eTCP,StopFlag,ActiveFlag);
                            if StopFlag == 1 || GUI_Variables.Stopped == 1
                                GUI_Variables.Stopped = 0;
                                StopFlag = 0;
                                fwrite(eTCP,'done');
                                set(handles.StatusText,'String',...
                                    horzcat(['Optimization stopped by user.'...
                                    ' To restart mid-generation select checkbox'...
                                    ' and fill last completed condition.'...
                                    ' Zero torque control is active. Last condition'...
                                    ' completed is Condition ',num2str(ConditionNumber-1)]));
                                if NumberofParams == 1
                                    set(handles.TRQ_SetpointText,'String',...
                                         num2str(0));
                                    SetpointHistory = [SetpointHistory; [0]];
                                    set(handles.TRQSetpointHistoryText,'String',...
                                         num2str(SetpointHistory(:,1)));
                                elseif NumberofParams == 2
                                    set(handles.TRQ_SetpointText,'String',...
                                         num2str(0));
                                    set(handles.RT_SetpointText,'String',...
                                         num2str(0));
                                    SetpointHistory = [SetpointHistory; [0 0]];
                                    set(handles.TRQSetpointHistoryText,'String',...
                                         num2str(SetpointHistory(:,1)));
                                    set(handles.RTSetpointHistoryText,'String',...
                                         num2str(SetpointHistory(:,2)));
                                end
                                GUI_Variables.TestDate = date;
                                error('User ended optimization.');
                            end
                            if ConditionNumber <= ConditionsPerGen && GUI_Variables.SeedNextGen == 0
                                disp('Controller should change to next set of parameters. If you wish to stop mid-generation, press ctrl c')
                                i = 0; %Initialize the counter for storing breaths
                                conditiondone = 0; %Initialize to make sure it knows condition isn't done yet.
                                ParamsForCondition = Params_full(ConditionNumber, :); %First row of params_full
                                fullrate = zeros(1,200);
                                breathtimes = zeros(1,200);

                            else %So if you are done with all the conditions for that generation
                                if GUI_Variables.SeedNextGen == 1 && ConditionNumber-1 < 3
                                    set(handles.StatusText,'String',...
                                        {'Need at least 3 conditions to seed next generation! Stopping optimization. To continue, start from mid-geneneration.'...
                                        ['Last completed condition was ',num2str(ConditionNumber-1),'.']} );
                                    GUI_Variables.SeedNextGen = 0;
                                    stop(TimerVar);
                                    error('Not enough conditions completed to seed next generation. Start from mid-generation');
                                else
                                
                                    orderedconds = ordering_conditions(SSdata); %(this is a conditionspergen rows by 2+params columns matrix.

                                    %Scale back into from 0 to 1, 0 to 1 for both parameters for use in
                                    %CMA. 
                                    %orderedconds(:,1) = orderedconds(:,1);
                                    %orderedconds(:,2) = orderedconds(:,2);
                                    if NumberofParams == 1
                                        orderedconds(:,3) = orderedconds(:,3)/Peak_torque;
                                        Next_Gen = @create_next_generation_newCMA_1Param;
                                    elseif NumberofParams == 2
                                        orderedconds(:,3) = orderedconds(:,3)/Peak_torque;
                                        orderedconds(:,4) = orderedconds(:,4)/100; 
                                        Next_Gen = @create_next_generation_newCMA_2Param;
                                    end
                                    Old_Params_full = Params_full; %Save old parameters for plotting

                                    [xmean, mu, weights, ps, pc, C, c1, cmu, cc, sigma, cs, damps,...
                                        chiN, eigeneval, invsqrtC, counteval, x, lambda, Params_full, N, B, D, mueff] = ...
                                        Next_Gen(orderedconds, xmean, mu, weights,ps, pc, C, c1,...
                                        cmu, cc, sigma, cs, damps, chiN, eigeneval, invsqrtC, counteval, x, lambda, N, B, D, mueff, Peak_torque, Min_torque)

                                    %When params_full comes out of the CMA, it comes out as 0 to 1, 0
                                    %to 1. 
                                    %Turn Params_full into actual numbers.
                                    if NumberofParams == 1
                                        Params_full(:,1) = Params_full(:,1)*Peak_torque; %Peak Torque (Nm)
                                        xmean_num = xmean'.* [Peak_torque];
                                    elseif NumberofParams == 2
                                        Params_full(:,1) = Params_full(:,1)*Peak_torque; %Peak Torque (Nm)
                                        Params_full(:,2) = Params_full(:,2)*100;  %Rise time (% of stance time)
                                        xmean_num = xmean'.* [Peak_torque, 100];
                                    end
                                    generationdone=1;
                                    stop(TimerVar) %Stops the timer. 
                                end

                            end
                        end
                        break
                    else
                    end
                end
            elseif generationdone == 1
                break
            elseif ActiveFlag == 0
            end
        else
            break
        end
    end

    disp('The generation is complete.')
    disp('The next set of parameters is')
    Params_full
    disp('The new xmean is')
    xmean_num
    xmean_history = [xmean_history; xmean_num]
    set(handles.StatusText,'String',...
        {'The generation is complete.',...
        'The next set of parameters is',...
        num2str(Params_full),...
        'The new xmean is',...
        num2str(xmean_num)});
    if NumberofParams == 1
        set(handles.XMeanTrq,'String',...
            num2str(xmean_history(:,1)));
    elseif NumberofParams == 2
        set(handles.XMeanTrq,'String',...
            num2str(xmean_history(:,1)));
        set(handles.XMeanRT,'String',...
           num2str(xmean_history(:,2)));
    end

    %Save all of the variables so you can load them in for the next generation.

    clearvars GenerationAcc 
    VarsToIgnore = ['eventdata|GUI_Variables|handles|hObject|ActiveFlag'...
        '|DoneFirstValue|SendValueFlag|SetpointHistory|eTCP'];
    save(fullfile(saveDir,['Completion_of_Gen_', num2str(GenerationNumber),'_',SSID]),...
        '-regexp',['^(?!',VarsToIgnore,'$).']);
    set(handles.StatusText,'String',...
        ['Saved file "Completion_of_Gen_',num2str(GenerationNumber),'_',SSID,'"']);

    fwrite(eTCP,'done') %Send completion code to exo GUI

    set(handles.StatusText,'String',...
        ['Concluded optimizer generation. To start next generation adjust',...
        ' optimizer settings and press "Start Optimizer"']);
    set(handles.GenNumber,'String',num2str(GenerationNumber+1));
    set(handles.GenNumber,'Value',GenerationNumber+1);
    set(handles.MidGenCheckbox,'Value',0);
    set(handles.LastConditionCompleted,'enable','off');
    set(handles.LastConditionCompleted,'String','0');
    set(handles.TestDate,'enable','off');
    set(handles.TestDate,'String','DD-Mon-YYYY');
    GUI_Variables.GenNum = get(handles.GenNumber,'Value');
    GUI_Variables.MidGen = get(handles.MidGenCheckbox,'Value');
    GUI_Variables.CompleteCond = str2double(get(handles.LastConditionCompleted,'String'));
    GUI_Variables.TestDate = date;
    GUI_Variables.pullDir = [GUI_Variables.SSID,'_',GUI_Variables.TestDate];
    GUI_Variables.SeedNextGen = 0;
    GUI_Variables.Streaming = 0;
    
    currentDir = cd;                                % Current directory
    saveDir = [currentDir,'\',SSID,'_',date,'\'];   % Save directory
    
    figure(100);
    hold on
    plot(1:GenerationNumber,xmean_history(:,1)/Subject_mass,'-k','LineWidth',2);
    plot(1:GenerationNumber,xmean_history(:,1)/Subject_mass,'xr','MarkerSize',10);
    GenPlot = repmat(1:GenerationNumber,length(Old_Params_full),1);
    for j = 1:length(Old_Params_full)
        plot(GenPlot(:,GenerationNumber),Old_Params_full/Subject_mass,'*k','MarkerSize',5);
    end
    xlabel('Generation Number')
    ylabel('Normalized Torque Magnitude [N*m/kg]')
    title('Mean Torque Setpoint Convergence')
    set(gcf,'Name',horzcat('Trq_Convergence_',SSID))
    set(gca,'XTick',1:GenerationNumber);
    set(gca,'FontName','Arial');
    set(gca,'FontWeight','Bold');
    hold off
    
    saveas(gcf,fullfile(saveDir,[get(gcf,'Name'),'.png']));
    
    if NumberofParams==2
        figure(101);
        plot(1:GenerationNumber,xmean_history(:,2),'-k','LineWidth',2);
        hold on
        plot(1:GenerationNumber,xmean_history(:,2),'xr','MarkerSize',10);
        xlabel('Generation Number')
        ylabel('Rise Time Percentage [% Stance]')
        title('Rise Time Setpoint Convergence')
        set(gcf,'Name',horzcat('RiseTime_Convergence_',SSID))
        set(gca,'XTick',1:GenerationNumber);
        set(gca,'FontName','Arial');
        set(gca,'FontWeight','Bold');
        hold off

       saveas(gcf,fullfile(saveDir,[get(gcf,'Name'),'.png']));
    end
    
    close(figure(102));
    figure(102);
    LineColors = {'-b','-g','-r','-m','-k','-c','-y',...
        '--b','--g','--r','--m','--k','--c','--y'};
    PointColors = {'xb','xg','xr','xm','xk','xc','xy',...
        'ob','og','or','om','ok','oc','oy'};
    FullMet=Full_Metabolic_Data_to_Save;
    hold on;
    for i = 1:length(FullMet)
        plot(FullMet{1,i}(:,2),Full_y_bar{1,i}',LineColors{i}); %Plot fit metabolic data
        LegendStr{i} = [num2str(Old_Params_full(i)),' Nm'];
    end
    for i = 1:length(FullMet)
        plot(FullMet{1,i}(:,2),FullMet{1,i}(:,1),PointColors{i}); %Plot raw metabolic data
    end
    legend(LegendStr,'Location','Best');
    xlabel('Condition Time [s]')
    ylabel('Metabolic Rate [W]')
    title(['Metabolic Rate Comparison for ','Gen ',num2str(GenerationNumber)]);
    set(gcf,'Name',['MetRateComparison_',SSID,'_Gen_',num2str(GenerationNumber)]);
    set(gca,'FontName','Arial');
    set(gca,'FontWeight','Bold');
    hold off
    
    saveas(gcf,fullfile(saveDir,[get(gcf,'Name'),'.png']));      
        
end    
    
end

% --- Timer function callback
function updater(obj, event)

global GUI_Variables

A = get(obj, 'UserData');
A(1, :) = A(1, :)+1; %incrementing i
i = A(1,1);
if i+2>length(A(:,1)) %Because that's all your remaining conditions
   %Option 1
    if GUI_Variables.NumParams == 1
        A(2,:) = [0]; %Set the controller to 0, or at least set this variable 
    %to 0 and then a line can be added to exoskeleton code to go to 
    %zero impedance mode if this happens.
    %Option 2
    %A(2,:) = A(2,:); %Alterantively, you can leave it on the last
    %controller while it calculates the next generation before you turn off
    %the assistance and have the participant stop walking. 
    elseif GUI_Variables.NumParams == 2
        A(2,:) = [0,0];
    end
else
    A(2, :) = A(i+2, :); %Updates row 2 to be the controller you should use. 
end
set(obj, 'UserData', A)
end

% --- Check Value Function
% Checks for values written from the exo computer to start/pause the
% optimization
function [ValueFlag,StopFlag,ActiveFlag] = CheckValue(tTCP,StopFlag,ActiveFlag)
%Check the value stored in the buffer
ValueFlag = 1;
if tTCP.BytesAvailable>0 && StopFlag ~= 1
    str = fread(tTCP,tTCP.BytesAvailable);
    if char(str') == "start"
        ValueFlag = 1;
        StopFlag = 0;
        ActiveFlag = 1;
    elseif char(str') == "end"
        ValueFlag = 0;
        StopFlag = 1;
        ActiveFlag = 0;
        return
    else
        ValueFlag = 1;
        StopFlag = 0;
    end
elseif tTCP.BytesAvailable == 0
    return
end

end
    
% --- Executes on button press in StopOptimizer.
function StopOptimizer_Callback(hObject, eventdata, handles)
% hObject    handle to StopOptimizer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

GUI_Variables.Stopped = 1;

set(handles.StatusText,'String',...
    'Stopping optimization and saving last completed condition.');
disp('Program terminating...');
end


function SubjectMass_Callback(hObject, eventdata, handles)
% hObject    handle to SubjectMass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

try
    GUI_Variables.SubjectMass = str2double(get(hObject,'String'));
    if isnan(GUI_Variables.SubjectMass)
        error();
    end
    set(handles.StatusText,'String',...
    horzcat('Subject mass entered: ',...
    num2str(GUI_Variables.SubjectMass),' kg'));    
catch
    set(handles.StatusText,'String',...
        'Invalid mass input. Try again.');
end
end

% Hints: get(hObject,'String') returns contents of SubjectMass as text
%        str2double(get(hObject,'String')) returns contents of SubjectMass as a double


% --- Executes during object creation, after setting all properties.
function SubjectMass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SubjectMass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function MaxPeakTorque_Callback(hObject, eventdata, handles)
% hObject    handle to MaxPeakTorque (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

try
    GUI_Variables.PkTRQ = str2double(get(hObject,'String'));
    if isnan(GUI_Variables.PkTRQ)
        error();
    end
    set(handles.StatusText,'String',...
    horzcat('Peak torque entered: ',...
    num2str(GUI_Variables.PkTRQ),' N*m/kg'));
catch
    set(handles.StatusText,'String',...
        'Invalid peak torque input. Try again');
end
% Hints: get(hObject,'String') returns contents of MaxPeakTorque as text
%        str2double(get(hObject,'String')) returns contents of MaxPeakTorque as a double
end


% --- Executes during object creation, after setting all properties.
function MaxPeakTorque_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxPeakTorque (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function SubjectID_Callback(hObject, eventdata, handles)
% hObject    handle to SubjectID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

try
    GUI_Variables.SSID = get(hObject,'String');
    set(handles.StatusText,'String',...
    horzcat('Subject ID entered: ',GUI_Variables.SSID)); 
catch
    set(handles.StatusText,'String',...
        'Invalid subject ID input. Try again.');
end

if ~any(isnan(GUI_Variables.SSID)) && GUI_Variables.GenNum>1
    set(handles.TestDate,'enable','on')
end
% Hints: get(hObject,'String') returns contents of SubjectID as text
%        str2double(get(hObject,'String')) returns contents of SubjectID as a double
end


% --- Executes during object creation, after setting all properties.
function SubjectID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SubjectID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function GenNumber_Callback(hObject, eventdata, handles)
% hObject    handle to GenNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GenNumber as text
%        str2double(get(hObject,'String')) returns contents of GenNumber as a double

global GUI_Variables

try
    GUI_Variables.GenNum = abs(str2double(get(hObject,'String')));
    if isnan(GUI_Variables.GenNum)
        error();
    end
    set(handles.StatusText,'String',...
    horzcat('Set generation number: ',...
    num2str(GUI_Variables.GenNum)));
catch
    set(handles.StatusText,'String',...
        'Invalid generation number. Try again.');
end

if ~any(isnan(GUI_Variables.SSID)) && (GUI_Variables.GenNum>1 || GUI_Variables.MidGen==1)
    set(handles.TestDate,'enable','on')
elseif GUI_Variables.GenNum == 1 && GUI_Variables.MidGen == 0
    set(handles.TestDate,'enable','off')
    set(handles.TestDate,'String','DD-Mon-YYYY')
    GUI_Variables.pullDir = ' ';
    GUI_Variables.TestDate = ' ';
end

end




% --- Executes during object creation, after setting all properties.
function GenNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GenNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function LastConditionCompleted_Callback(hObject, eventdata, handles)
% hObject    handle to LastConditionCompleted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

try
    if GUI_Variables.MidGen == 1
        GUI_Variables.CompleteCond = abs(str2double(get(hObject,'String')));
        if isnan(GUI_Variables.CompleteCond)
            set(handles.StatusText,'String',...
                'Invalid last completed condition entered. Try again.');
        end
        set(handles.StatusText,'String',...
        horzcat('Last condition completed: ',...
        num2str(GUI_Variables.CompleteCond)));
    else
        disp('Select mid generation textbox to define last completed condition.');
        set(handles.StatusText,'String',...
            'Select mid generation textbox to define last completed condition.');
    end
catch
    set(handles.StatusText,'String',...
        'Invalid number of completed conditions. Try again');
end
end

% Hints: get(hObject,'String') returns contents of LastConditionCompleted as text
%        str2double(get(hObject,'String')) returns contents of LastConditionCompleted as a double


% --- Executes during object creation, after setting all properties.
function LastConditionCompleted_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LastConditionCompleted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created untl after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function ConditionTime_Callback(hObject, eventdata, handles)
% hObject    handle to ConditionTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

try
    GUI_Variables.ConditionTime = str2double(get(hObject,'String'));
    set(handles.ConditionTime,'String',num2str(GUI_Variables.ConditionTime));
    set(handles.StatusText,'String',...
        horzcat('Input time between conditions: ',...
        num2str(GUI_Variables.ConditionTime),...
        ' seconds'));
    if isnan(GUI_Variables.ConditionTime)
        error();
    end
catch
    set(handles.StatusText,'String',...
        'Invalid time between conditions entered. Try again and ensure input in seconds.');
end
end

% Hints: get(hObject,'String') returns contents of ConditionTime as text
%        str2double(get(hObject,'String')) returns contents of ConditionTime as a double


% --- Executes during object creation, after setting all properties.
function ConditionTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ConditionTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in OneParamCheck.
function OneParamCheck_Callback(hObject, eventdata, handles)
% hObject    handle to OneParamCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OneParamCheck

global GUI_Variables

if get(hObject,'Value')
    GUI_Variables.NumParams = 1;
    GUI_Variables.CondPerGen = 5;
    set(handles.TwoParamCheck,'Value',0);
end

end


% --- Executes on button press in TwoParamCheck.
function TwoParamCheck_Callback(hObject, eventdata, handles)
% hObject    handle to TwoParamCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TwoParamCheck

global GUI_Variables

if get(hObject,'Value')
    GUI_Variables.NumParams = 2;
    GUI_Variables.CondPerGen = 8;
    set(handles.OneParamCheck,'Value',0);
end

end

function MinTRQ_Callback(hObject, eventdata, handles)
% hObject    handle to MinTRQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinTRQ as text
%        str2double(get(hObject,'String')) returns contents of MinTRQ as a double

global GUI_Variables

try 
    GUI_Variables.MinTRQ = str2double(get(hObject,'String'));
    set(handles.StatusText,'String',...
        ['Min torque entered: ', get(hObject,'String'), ' N*m/kg']);
catch
    set(handles.StatusText,'String',...
        'Invalid minimum torque defined. Try again');
end
    
end

% --- Executes during object creation, after setting all properties.
function MinTRQ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinTRQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function TestDate_Callback(hObject, eventdata, handles)
% hObject    handle to TestDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TestDate as text
%        str2double(get(hObject,'String')) returns contents of TestDate as a double

% --- Executes during object creation, after setting all properties.

global GUI_Variables

try
    GUI_Variables.TestDate = get(hObject,'String');
    set(handles.TestDate,'String',GUI_Variables.TestDate);
    set(handles.StatusText,'String',...
        ['Input last test date: ', GUI_Variables.TestDate]);
    Months = {'Jan','Feb','Mar','Apr','May','Jun','Jul',...
        'Aug','Sep','Oct','Nov','Dec'}; %Valid months
    DirFiles = dir; %Get all directory files
    ValidDir = [];
    for i = 1:length(DirFiles)
        if DirFiles(i).isdir && contains(DirFiles(i).name,GUI_Variables.SSID)
            ValidDir = [ValidDir, i];
        else
        end
    end
    for i = 1:length(ValidDir)
        CheckDir{i} = DirFiles(ValidDir(i)).name; %Get all directories containing subject ID
    end
    CheckDate = strsplit(GUI_Variables.TestDate,'-');
    if length(CheckDate) ~= 3 %If not three distinct cells throw error
        set(handles.StatusText,'String',...
            ['Invalid date format. Need DD-Mon-YYYY']);
        error();
    elseif length(CheckDate{1}) ~= 2 %If not two days indicated throw error
        set(handles.StatusText,'String',...
           ['Invalid date format. Need DD']);
        error();
    elseif str2double(CheckDate{1}) < 1 || str2double(CheckDate{1}) > 31 %if not valid day
        set(handles.StatusText,'String',...
            ['Invalid date format. Day definition outside of typical month']);
        error();
    elseif ~any(strcmp(Months,CheckDate{2})) %If not valid month
        set(handles.StatusText,'String',...
            ['Invalid date format. Check month definition. Should be "Mon"']);
        error();
    elseif length(CheckDate{3}) ~= 4 %If not four number year defined
        set(handles.StatusText,'String',...
            ['Invalid date format. Check year definition. Should be YYYY']);
        error();        
    elseif str2double(CheckDate{3}) < 0 %If negative year defined
        set(handles.StatusText,'String',...
            ['Invalid date format. Check year']);
        error();
    elseif ~any(contains(CheckDir,GUI_Variables.TestDate))
        set(handles.StatusText,'String',...
            ['Chosen directory not found! Check SSID and date']);
        error();
    elseif strcmp(GUI_Variables.TestDate,date) %If current date
        set(handles.StatusText,'String',...
            ['Entered date is current date! No need to re-define.']);
        set(handles.TestDate,'String','DD-Mon-YYYY');
        GUI_Variables.TestDate = ' ';
    else
        GUI_Variables.pullDir = [GUI_Variables.SSID,'_',GUI_Variables.TestDate];
        set(handles.StatusText,'String',...
            ['Defined directory to pull from: ',GUI_Variables.pullDir]); 
    end
catch
end

end
    
function TestDate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TestDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in SeedNextGen.
function SeedNextGen_Callback(hObject, eventdata, handles)
% hObject    handle to SeedNextGen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GUI_Variables

if GUI_Variables.Streaming == 1 %If we're actively sending data, seed after completion of current condition

    GUI_Variables.SeedNextGen = 1;
    set(handles.StatusText,'String','Attempting to seed next generation after completing current condition...');

else %If we're stopped we can attempt to seed from the available saved conditions
    
    GenerationNumber = GUI_Variables.GenNum;
    NextConditionMidGen = GUI_Variables.CompleteCond;
    SSID = GUI_Variables.SSID;
    
    currentDir = cd;                                % Current directory
    saveDir = [currentDir,'\',SSID,'_',date,'\'];   % Save directory
    mkdir(saveDir);
    
    set(handles.StatusText,'String','Attempting to seed next generation from saved condition...');
    if strcmp(GUI_Variables.TestDate,date) || strcmp(GUI_Variables.TestDate,' ')  %Compare current date to entered date
        load(fullfile(saveDir,['All_Saved_Data_Following_Gen_', num2str(GenerationNumber), '_Cond_', num2str(NextConditionMidGen),'_',SSID]))
    else
        load(fullfile(GUI_Variables.pullDir,['All_Saved_Data_Following_Gen_', num2str(GenerationNumber), '_Cond_', num2str(NextConditionMidGen),'_',SSID]))
    end
    currentDir = cd;                                % Current directory
    saveDir = [currentDir,'\',SSID,'_',date,'\'];   % Save directory
    
    orderedconds = ordering_conditions(SSdata); %(this is a conditionspergen rows by 2+params columns matrix.

    %Scale back into from 0 to 1, 0 to 1 for both parameters for use in
    %CMA. 
    %orderedconds(:,1) = orderedconds(:,1);
    %orderedconds(:,2) = orderedconds(:,2);
    if NumberofParams == 1
        orderedconds(:,3) = orderedconds(:,3)/Peak_torque;
        Next_Gen = @create_next_generation_newCMA_1Param;
    elseif NumberofParams == 2
        orderedconds(:,3) = orderedconds(:,3)/Peak_torque;
        orderedconds(:,4) = orderedconds(:,4)/100; 
        Next_Gen = @create_next_generation_newCMA_2Param;
    end
    Old_Params_full = Params_full; %Save old parameters for plotting

    [xmean, mu, weights, ps, pc, C, c1, cmu, cc, sigma, cs, damps,...
        chiN, eigeneval, invsqrtC, counteval, x, lambda, Params_full, N, B, D, mueff] = ...
        Next_Gen(orderedconds, xmean, mu, weights,ps, pc, C, c1,...
        cmu, cc, sigma, cs, damps, chiN, eigeneval, invsqrtC, counteval, x, lambda, N, B, D, mueff, Peak_torque, Min_torque)
    
    if NumberofParams == 1
        Params_full(:,1) = Params_full(:,1)*Peak_torque; %Peak Torque (Nm)
        xmean_num = xmean'.* [Peak_torque];
    elseif NumberofParams == 2
        Params_full(:,1) = Params_full(:,1)*Peak_torque; %Peak Torque (Nm)
        Params_full(:,2) = Params_full(:,2)*100;  %Rise time (% of stance time)
        xmean_num = xmean'.* [Peak_torque, 100];
    end
    
    disp('The generation is complete.')
    disp('The next set of parameters is')
    Params_full
    disp('The new xmean is')
    xmean_num
    xmean_history = [xmean_history; xmean_num]
    set(handles.StatusText,'String',...
        {'The generation is complete.',...
        'The next set of parameters is',...
        num2str(Params_full),...
        'The new xmean is',...
        num2str(xmean_num)});
    if NumberofParams == 1
        set(handles.XMeanTrq,'String',...
            num2str(xmean_history(:,1)));
    elseif NumberofParams == 2
        set(handles.XMeanTrq,'String',...
            num2str(xmean_history(:,1)));
        set(handles.XMeanRT,'String',...
           num2str(xmean_history(:,2)));
    end
    
    clearvars GenerationAcc 
    VarsToIgnore = ['eventdata|GUI_Variables|handles|hObject|ActiveFlag'...
        '|DoneFirstValue|SendValueFlag|SetpointHistory|eTCP'];
    
    save(fullfile(saveDir,['Completion_of_Gen_', num2str(GenerationNumber),'_',SSID]),...
        '-regexp',['^(?!',VarsToIgnore,'$).']);
    set(handles.StatusText,'String',...
        ['Saved file "Completion_of_Gen_',num2str(GenerationNumber),'_',SSID,'"']);
    
    set(handles.StatusText,'String',...
        ['Concluded optimizer generation. To start next generation adjust',...
        ' optimizer settings and press "Start Optimizer"']);
    set(handles.GenNumber,'String',num2str(GenerationNumber+1));
    set(handles.GenNumber,'Value',GenerationNumber+1);
    set(handles.MidGenCheckbox,'Value',0);
    set(handles.LastConditionCompleted,'enable','off');
    set(handles.LastConditionCompleted,'String','0');
    set(handles.TestDate,'enable','off');
    set(handles.TestDate,'String','DD-Mon-YYYY');
    GUI_Variables.GenNum = get(handles.GenNumber,'Value');
    GUI_Variables.MidGen = get(handles.MidGenCheckbox,'Value');
    GUI_Variables.CompleteCond = str2double(get(handles.LastConditionCompleted,'String'));
    GUI_Variables.TestDate = date;
    GUI_Variables.pullDir = [GUI_Variables.SSID,'_',GUI_Variables.TestDate];
    GUI_Variables.SeedNextGen = 0;
    GUI_Variables.Streaming = 0;
    
    figure(100);
    hold on
    plot(1:GenerationNumber,xmean_history(:,1)/Subject_mass,'-k','LineWidth',2);
    plot(1:GenerationNumber,xmean_history(:,1)/Subject_mass,'xr','MarkerSize',10);
    GenPlot = repmat(1:GenerationNumber,length(Old_Params_full),1);
    for j = 1:length(Old_Params_full)
        plot(GenPlot(:,GenerationNumber),Old_Params_full/Subject_mass,'*k','MarkerSize',5);
    end
    xlabel('Generation Number')
    ylabel('Normalized Torque Magnitude [N*m/kg]')
    title('Mean Torque Setpoint Convergence')
    set(gcf,'Name',horzcat('Trq_Convergence_',SSID))
    set(gca,'XTick',1:GenerationNumber);
    set(gca,'FontName','Arial');
    set(gca,'FontWeight','Bold');
    hold off
    
    saveas(gcf,fullfile(saveDir,[get(gcf,'Name'),'.png']));
    
    if NumberofParams==2
        figure(101);
        plot(1:GenerationNumber,xmean_history(:,2),'-k','LineWidth',2);
        hold on
        plot(1:GenerationNumber,xmean_history(:,2),'xr','MarkerSize',10);
        xlabel('Generation Number')
        ylabel('Rise Time Percentage [% Stance]')
        title('Rise Time Setpoint Convergence')
        set(gcf,'Name',horzcat('RiseTime_Convergence_',SSID))
        set(gca,'XTick',1:GenerationNumber);
        set(gca,'FontName','Arial');
        set(gca,'FontWeight','Bold');
        hold off

       saveas(gcf,fullfile(saveDir,[get(gcf,'Name'),'.png']));
    end
    
    close(figure(102));
    figure(102);
    LineColors = {'-b','-g','-r','-m','-k','-c','-y',...
        '--b','--g','--r','--m','--k','--c','--y'};
    PointColors = {'xb','xg','xr','xm','xk','xc','xy',...
        'ob','og','or','om','ok','oc','oy'};
    FullMet=Full_Metabolic_Data_to_Save;
    hold on;
    for i = 1:length(FullMet)
        plot(FullMet{1,i}(:,2),Full_y_bar{1,i}',LineColors{i}); %Plot fit metabolic data
        LegendStr{i} = [num2str(Old_Params_full(i)),' Nm'];
    end
    for i = 1:length(FullMet)
        plot(FullMet{1,i}(:,2),FullMet{1,i}(:,1),PointColors{i}); %Plot raw metabolic data
    end
    legend(LegendStr,'Location','Best');
    xlabel('Condition Time [s]')
    ylabel('Metabolic Rate [W]')
    title(['Metabolic Rate Comparison for ','Gen ',num2str(GenerationNumber)]);
    set(gcf,'Name',['MetRateComparison_',SSID,'_Gen_',num2str(GenerationNumber)]);
    set(gca,'FontName','Arial');
    set(gca,'FontWeight','Bold');
    hold off
    
    saveas(gcf,fullfile(saveDir,[get(gcf,'Name'),'.png']));      
    
end

end
