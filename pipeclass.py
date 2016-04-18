#coding=utf-8;
from math import pi;import math;
import numpy;
import scipy;
from scipy.optimize import fsolve
from iapws import IAPWS95;
from iapws import IAPWS97;

class pipe():
    #parameter(design)
    le=0.15;lc=1.15;la=0.1;leff=(le+lc)/2+la;l=le+lc+la;
    dv=0.121;
    di=0.151;
    do= 0.224;
    d=6.25e-5;
    N=7870.0;
    kp=300.0;kw=400.0
    rhoPipeWall=4200.0;rhoNet=3100.0
    cPipWall=46.0;cNet=46.0
    Tc=400.0;Te=460.0
    P=5.5;
    W=0.0001;
    theta=pi/4;

    Q=12.0
    fasshas=False
    def calStablPipeT(self):
        self.k=1/(1/(2*pi)*(1/self.kp*math.log(self.do/self.di)+1/self.kw*math.log(self.di/self.dv)))
        self.T=(self.Te*self.le+self.Tc*self.lc)/(self.lc+self.le)
        self.Q=(self.Te-self.T)/(1/self.k)


    def calStablPipeT2(self):
        self.k=1/(1/(2*pi)*(1/self.kp*math.log(self.do/self.di)+1/self.kw*math.log(self.di/self.dv)))
        self.T=self.Te-self.Q/self.le/self.k;
        print self.k,self.T,"************"


    '''水及水蒸气的参数计算'''
    def getRefrigProp(self):

        print self.T,"T"
        if self.T<274.0:
            print "冷冻"
            self.T=440.0
        if self.T>540.0:
            print "超临界水"
            self.T=440.0
        steam=IAPWS95(T=self.T,x=1);
        water=IAPWS95(T=self.T,x=0);
        #parameter in process
        self.hfg=water.Hvap*1000.0;
        self.sigma=water.sigma
        self.Pv=steam.P*10e6;
        self.uv=steam.mu
        self.rhov=steam.rho;

        self.gammav=steam.cp0/steam.cv0;
        self.Rv=steam.s0
        self.rhol=water.rho
        self.Ul=water.nu*self.rhol
        self.Uv=steam.nu*self.rhov
        self.cpv=steam.cp0
        self.cpwater=water.cp0

    def calCriticlQ(self):
        #黏性极限的计算
        self.Av=pi*self.dv**2/4;
        self.Qvmax=self.dv**2*self.hfg/64/self.uv/self.leff*self.rhov*self.Av;
        #携带极限
        self.rns=1/(2.0*self.N)-self.d/2
        self.Qemax=(self.rhov*self.sigma/2/self.rns)**0.5*self.hfg*self.Av
        #声速极限
        self.Qsmax=self.Av*self.rhov*self.hfg*(self.gammav*self.Rv*self.T/2/(self.gammav+1))
        #毛细极限的计算
        self.rc=1/(self.N*2)
        self.epsilon=1-1.05*pi*self.N*self.d/4
        self.K=self.d**2*self.epsilon**2/122.0/(1-self.epsilon**2)
        self.Aw=pi*(self.di**2-self.dv**2)/4
        self.rhv=self.dv/2
        Fl=self.Ul/(self.K*self.Aw*self.rhol*self.hfg)
        Fv=16.0*self.Uv/(2*self.Av*self.rhv**2*self.rhov*self.hfg)
        Qcmax=2*self.sigma/(self.rc*(Fl+Fv)*self.leff)
        RevCri=2*self.rhv*Qcmax/self.Av/self.Uv/self.hfg;
        Mv=Qcmax/(self.Av*self.rhov*self.hfg*(self.gammav*self.Rv*self.T)**0.5)
        _der=self.Uv/(self.Av*self.rhv**2*self.rhov*self.hfg)
          #私有方法校验计算雷诺数和马赫数 书本54页(Rev<=2300?Mv<=0.2?)
        _mthod={(True,True):(lambda :8*_der),
                (True,False):(lambda :8*_der*(1+(self.gammav-1)/2*Mv**2)**(-0.5)),
                (False,True):(lambda :0.038*_der*(2*self.rhv*Qcmax/self.Av/self.hfg/self.Uv)**(3/4)),
                (False,False):(lambda :0.038*_der*(2*self.rhv*Qcmax/self.Av/self.hfg/self.Uv)**(3/4)*(1+(self.gammav-1)/2*Mv**2)**(-0.5))}
        while True:
            _ss=(RevCri<=2300,Mv<=0.2)
            Fv=_mthod[_ss]()
            Qcmax=2*self.sigma/(self.rc*(Fl+Fv)*self.leff)
            RevCri=2*self.rhv*Qcmax/self.Av/self.Uv/self.hfg;
            Mv=Qcmax/(self.Av*self.rhov*self.hfg*(self.gammav*self.Rv*self.T)**0.5)
            if _ss==(RevCri<=2300,Mv<=0.2):
                break
            else:
                pass
        self.Qcmax=Qcmax
        #计算热管质量
        mPipWall=pi*(self.do**2-self.di**2)/4*self.rhoPipeWall*self.l
        mSteam=pi*self.dv**2/4.0*self.rhov*self.l+pi*(self.di**2-self.dv**2)/4*(self.leff-self.la)*self.rhov*self.epsilon
        mWater=pi*(self.di**2-self.dv**2)/4*self.leff*self.rhol*self.epsilon
        mNet=pi*(self.di**2-self.dv**2)/4*self.l*self.rhoNet*(1-self.epsilon)
        self.mTotal=mNet+mSteam+mWater+mPipWall
        #冷启动极限

        _cThemalcpPerM=(mNet*self.cNet+mSteam*self.cpv+mWater*self.cpwater+mPipWall*self.cPipWall)/self.l
        self.Q1coldmax=self.epsilon*self.rhol*self.Aw*self.hfg/_cThemalcpPerM/(self.T-10)
        print self.Q1coldmax,"self.Q1coldmax"
    #压降线性方程
    def _calMmax(self,mmax):
        detaPmao=2*self.sigma*math.cos(self.theta)/self.rc
        detaPvsc=1/self.dv*(462.0*self.T*2*pi)**(0.5)*(mmax)/pi/self.lc+detaPmao
        detaPvse=1/self.dv*(462.0*self.T*2*pi)**(0.5)*(mmax)/pi/self.le
        #下面阻力
        Rev=mmax*self.dv/self.Uv/(pi*self.dv**2/4)
        detaPl=self.Ul*mmax*self.leff/self.rhol/self.K/self.Aw/self.epsilon
        if Rev<1:
            detaPv=4*self.Uv*self.l*mmax*self.hfg/pi/self.rhov/self.rhv**4/self.hfg
        else:
            detaPv=(1-4/pi**2)*(mmax*self.hfg)**2/8/self.rhov/self.rhv**4/self.hfg**2
        #print self.Uv,Rev,self.la,self.rhov,self.rhv
        detaPla=0.0655*self.Uv**2.0*Rev**(1.75)*self.la/self.rhov/(self.rhv**3)
        y=6.0
        #print y+detaPl;
        #print y.__class__,mmax,"detaPmao",detaPmao.__class__,'detaPvsc',detaPvsc-detaPmao,'detaPvse',detaPvse.__class__,'detaPl',detaPl.__class__,'detaPv',detaPv.__class__,'detaPla',detaPla.__class__,'Rev',Rev.__class__
        zero=detaPvsc+detaPvse-detaPl-detaPv-detaPla+1.0-1.0
        return zero

    def calPipeWokingState(self):
        self.mmax=fsolve(self._calMmax,0.00000001)

    def calR(self):
        self.R2=math.log(self.do/self.di)/2.0/pi/self.kp/self.le
        self.R8=math.log(self.do/self.di)/2.0/pi/self.kp/self.lc
        self.R3=math.log(self.di/self.dv)/2.0/pi/self.kw/self.le
        self.R7=math.log(self.di/self.dv)/2.0/pi/self.kw/self.lc
        self.R10=4*self.leff/(self.kp*pi*(self.do**2-self.di**2))
        R11=4*self.leff/(self.kw*pi*(self.di**2-self.dv**2))
        self.R4=462.0*self.T*(2*pi*462.0*self.T)**(1/2)/self.l**2/self.Pv/pi/self.dv/self.le
        self.R6=462.0*self.T*(2*pi*462.0*self.T)**(1/2)/self.l**2/self.Pv/pi/self.dv/self.lc
        R5=128.0*self.Uv*self.leff*self.T/pi/self.rhov**2*self.dv**4*self.l**2
        return self.bingliandianzu(self.bingliandianzu(self.R3+self.R4+R5+self.R6+self.R7,R11)+self.R2+self.R8,self.R10)


    def bingliandianzu(self,r1,r2):
        return 1.0/(1.0/r1+1.0/r2);

    def calPipemmax(self):
        def _fun(Tcandmmax):
            Tc=Tcandmmax[0]
            mmax=Tcandmmax[1]


            self.T=(self.Te-Tc-(mmax*self.hfg-(self.Te-Tc)/self.R10)*(self.R2+self.R8))/(1+(self.R4+self.R3)/(self.R6+self.R7))+Tc+(mmax*self.hfg-(self.Te-Tc)/self.R10)*self.R8
            self.reset()
            self.calR()
            #print self.T,"55555555555544444444444444"
            mmax=Tcandmmax[1]

            #print self.Te-mmax*self.hfg/self.calR()-Tc,"??????",self._calMmax(mmax)
            print self.calR(),mmax*self.hfg,Tc

            return [self.Te-mmax*self.hfg*self.calR()-Tc,self._calMmax(mmax)+1.0-1.0]
        Tcandmmax=fsolve(_fun,[self.Tc,(self.Te-self.Tc)/self.calR()/self.hfg])
        self.mmax=Tcandmmax[1]
        self.Tc=Tcandmmax[0]
        #self.T=(self.Te-self.Tc-(self.mmax*self.hfg-(self.Te-self.Tc)/self.R10)*(self.R2+self.R8))/(1+(self.R4+self.R3)/(self.R6+self.R7))+self.Tc+(self.mmax*self.hfg-(self.Te-self.Tc)/self.R10)*self.R8
        #self.reset()
        self.Rrelative=self.calR()

    def init(self):
        self.calStablPipeT2()
        self.getRefrigProp()
        self.calCriticlQ()
        #self._calMmax()
        #self.calPipeWokingState()

    def reset(self):
        #self.calStablPipeT2()
        self.getRefrigProp()
        self.calCriticlQ()
        #self._calMmax()
        #self.calPipeWokingState()

pipe1=pipe()
pipe1.calStablPipeT2()
pipe1.reset()
pipe1.calPipemmax()
print pipe1.mmax,pipe1.hfg , pipe1.Q , pipe1.Qemax,"watching"
print "在以下设计处参数条件下：","蒸发段长",pipe1.le,"冷凝段长",pipe1.lc,"绝热段长",pipe1.la,"气腔直径",pipe1.dv ,"热管外径",pipe1.do,"热管内经",pipe1.di ,"毛细丝网直径",pipe1.d,"单位长度毛细数",pipe1.N   ,"热管外壁及吸液芯导热系数",pipe1.kp,pipe1.kw,"蒸发段温度",pipe1.Te ,"毛细材料接触角",pipe1.theta
print "热管的","Qcmax",pipe1.Qcmax,"Qvmax",pipe1.Qvmax,"Qemax",pipe1.Qemax,"Qsmax",pipe1.Qsmax,"热管冷凝段表面温度",pipe1.Tc,"热管内工质运行温度",pipe1.T,"热管蒸发段表面温度",pipe1.Te,"传热量",pipe1.mmax*pipe1.hfg,"等效热阻",pipe1.Rrelative,"冷启动极限是否发生",pipe1.Q1coldmax<1

