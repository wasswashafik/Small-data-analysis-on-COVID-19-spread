
function small_data_analysis

% The function makes data analysis of covid-19 in Kyrgyzstan, Kazakhstan
% and Russia. 
% SIR is based on Giovanni Valentini (2020). SIR Epidemic Spread Model
% (https://www.mathworks.com/matlabcentral/fileexchange/75100-sir-epidemic-spread-model).
% The data was collected from JHU
% https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series.
% SIR model is modified to allow beta to vary dynamically: the parameters
% are estimated to fit the data.



f = figure('Visible','off','Position',[360,500,450,285]);

hline   = uicontrol('Style','pushbutton',...
             'String','Linear fit','Position',[315,220,70,25],...
             'Callback',@linear,'FontSize',11);
hexp    = uicontrol('Style','pushbutton',...
             'String','Exponential fit','Position',[315,180,70,25],...
             'Callback',@exponential,'FontSize',11);
hsir = uicontrol('Style','pushbutton',...
             'String','SIR model','Position',[315,135,70,25],...
             'Callback',@sir,'FontSize',11);
htext  = uicontrol('Style','text','String','Select Data',...
           'Position',[325,90,60,15],'FontSize',11);
       
hpopup = uicontrol('Style','popupmenu',...
           'String',{'Kyrgyzstan','Kazakhstan','Russia'},...
           'Position',[300,50,100,25],...
           'Callback',@popup_menu_Callback,'FontSize',12);
ha = axes('Units','pixels','Position',[50,60,200,185]);
align([hline,hexp,hsir,htext,hpopup],'Right','None');


f.Units = 'normalized';
ha.Units = 'normalized';
hline.Units = 'normalized';
hexp.Units = 'normalized';
hsir.Units = 'normalized';
htext.Units = 'normalized';
hpopup.Units = 'normalized';

    
    data=load('data.mat');
    data=data.data;

    current_data=data.data_kg;
    
    fitdata=data.data_kg(198:end);
    plotdata=datetime(2020,10,1):datetime(2020,12,13);
    pl_n=datetime(2020,12,14):datetime(2020,12,20);
    Title='Kyrgyzstan';
    
    beta=data.beta_kg;
    N=6.3*10^6;
    T=427;
    plotdata_sir=datetime(2020,03,18):datetime(2021,05,18);
    plotdata_sir_obs=datetime(2020,03,18):datetime(2020,12,13);
    
f.Name = 'COVID-19 data analysis';

movegui(f,'center')

f.Visible = 'on';


   function popup_menu_Callback(source,eventdata) 
      
      str = get(source, 'String');
      val = get(source,'Value');
      
      switch str{val};
      case 'Kyrgyzstan' 
         current_data = data.data_kg;
         fitdata=data.data_kg(198:length(current_data));
         Title='Kyrgyzstan';
       
         beta=data.beta_kg;
         N=6.3*10^6;
          plotdata_sir=datetime(2020,03,18):datetime(2021,05,18);
    plotdata_sir_obs=datetime(2020,03,18):datetime(2020,12,13);
    T=427;
      case 'Kazakhstan' 
         current_data = data.data_kz;
         fitdata=data.data_kz(200:length(current_data));
         Title='Kazakhstan';
         
         
         beta=data.beta_kz;
         N=18.3*10^6;
          plotdata_sir=datetime(2020,03,16):datetime(2021,05,18);
    plotdata_sir_obs=datetime(2020,03,16):datetime(2020,12,13);
        T=429;
         
      case 'Russia' 
         current_data = data.data_ru;
         fitdata=data.data_ru(210:length(current_data));
         Title='Russia';
         beta=data.beta_ru;
         N=145*10^6;
         plotdata_sir=datetime(2020,03,06):datetime(2021,05,18);
        plotdata_sir_obs=datetime(2020,03,06):datetime(2020,12,13);
        T=439;
         
      end
   end



  function linear(source,eventdata) 
  
        nn=length(fitdata);
        obs=1:nn;
        linfit = fit(obs',fitdata','poly1');
        obs_n=nn+1:nn+7;
     
        
       cla;
       plot(plotdata, linfit(1:nn),'LineWidth',1.5);
       hold on;
       plot(pl_n,linfit(nn+1:nn+7),'--','Color','b');
       hold on;
       plot(plotdata, fitdata,'*','Color','black');
       hold off;
       legend('Linear fit','Projection','Observed','Location','NorthWest')
       title(Title);
       ylabel('confirmed cases');
  end

  function exponential(source,eventdata) 
 
        nn=length(fitdata);
        obs=1:nn;
        expfit = fit(obs',fitdata','exp1');
        obs_n=nn+1:nn+7;
    
        
       cla;
       plot(plotdata, expfit(1:nn),'LineWidth',1.5);
       hold on;
       plot(pl_n,expfit(nn+1:nn+7),'--','Color','b');
       hold on;
       plot(plotdata, fitdata,'*','Color','black');
       hold off;
       legend('Exponential fit','Projection','Observed','Location','NorthWest')
       title(Title);
        ylabel('confirmed cases');
  end

  function sir(source,eventdata) 
gamma=1/28;
delta = 1/90; 
I0 = 10; 
dt = 1; 

I = sir_model(beta(1,:),gamma,delta,N,I0,T,dt);
I2 = sir_model(beta(2,:),gamma,delta,N,I0,T,dt);
cla;

plot(plotdata_sir,I,'LineWidth',1.5);
xlim([datetime(2020,3,01) datetime(2021,05,30)]);
ylabel('confirmed cases');
hold on;

plot(plotdata_sir,I2,'LineWidth',1.5,'Color','green');

plot(plotdata_sir_obs,current_data,'.','Color','black');
       legend('SIR-current infection rate','SIR - slowed infection rate','Observed','Location','NorthWest')
       title(Title);
hold off;

  end

function I = sir_model(beta,gamma,delta,N,I0,T,dt)
  
    S = zeros(1,T/dt);
    S(1) = N;
    I = zeros(1,T/dt);
    I(1) = I0;
    R = zeros(1,T/dt);

    for tt = 1:(T/dt)-1
       
        dS = (-beta(tt)*I(tt)*S(tt) + delta*R(tt)) * dt;
        dI = (beta(tt)*I(tt)*S(tt) - gamma*I(tt)) * dt;
        dR = (gamma*I(tt) - delta*R(tt)) * dt;
        S(tt+1) = S(tt) + dS;
        I(tt+1) = I(tt) + dI;
        R(tt+1) = R(tt) + dR;
        
    end

end
end