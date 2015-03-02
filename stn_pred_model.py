from collections import namedtuple
from datetime import date as datetype
from datetime import time as timetype
from datetime import datetime, timedelta
from scipy.optimize import leastsq
from scipy.special import erf
from xml.dom import minidom
from xml.etree import ElementTree
import argparse
import ellipsoid
import math
import numpy as np
import os.path
import re
import sys

'''
stn_pred_model: 

A module for calculating functional representations of GNSS time series data, termed 
station coordinate prediction models. Each model consists of a set of terms such as 
constant velocity functions, step functions, cyclic functions etc. The model also 
holds information about when the time series has data, and about observations that 
have been excluded from the model calculation.

The module provides the following classes:

class model:  
    The main class encapsulating the station coordinate prediction model.
    This manages reading and saving the model as an XML file, read in time series data,
    fitting the model to the time series. 

class base_function: 
    The (abstract) base class for the function components of a model.  The following functions
    are derived from it.
    * offset
    * velocity
    * velocity_change
    * cyclic
    * annual
    * semiannual
    * equipment_offset
    * tectonic_offset
    * slow_slip
    * slow_slip_ramp
    * exponential_decay

class parameter:
    Base class for a parameter of a functional model, mainly concerned with conversion to
    and from the XML representation.  The following classes are derived from it:
    * offset_parameter
    * velocity_parameter
    * date_parameter

class exclude_obs:
    Class used to represent an observation excluded from a time series.

'''

refdate=datetime(2000,1,1)
daysperyear=365.25

cpm_tag='coordinate_prediction_model'
stn_tag='station'
outages_tag='outages'
outage_tag='outage'

                                    
dateformat='%d-%m-%Y'
datetimeformat='%Y-%m-%dT%H:%M:%S'

# Time conversion/utility functions

def asday( date ):
    '''
    Utility function to convert a time to a day number, referenced to refdate (2000-01-01)

    Input date can be 
    * a floating point number, treated as a day number
    * a string formatted as d-mm-yyyy hh:mm or yyyy-mm-ddThh:mm:ss
    * a python date or datetime class
    '''
    if type(date)==float:
        return date
    if type(date)==str:
        try:
            date=datetime.strptime(date,"%d-%m-%Y %H:%M")
        except:
            try:
                date=datetime.strptime(date,dateformat)
            except:
                date=datetime.strptime(date,datetimeformat)
    elif type(date) == datetype:
        date=datetime.combine(date,timetype(0,0,0))
    td = date.replace(tzinfo=None) - refdate
    return float(td.days)+float(td.seconds)/(60*60*24) 

def fromday( days ):
    '''
    Utility function to convert a day number (relative to refdate) to a python datetime object
    '''
    if type(days) == datetime:
        return days
    td = timedelta(days)
    return refdate+td

def days_array( dates ):
    '''
    Convert an array of dates to a numpy array of floating point day numbers. 
    
    Uses the asday function to convert the dates if necessary.
    '''
    if not isinstance(dates,np.ndarray):
        if not isinstance(dates,list):
            dates=[dates]
        dates=np.array(dates)
    if type(dates[0]) in (datetime, datetype):
        dates=np.array([asday(d) for d in dates])
    return dates


# Parameter types used to define model parameters

class parameter( object ):
    '''
    Base class for a parameter of a model function.

    Also useable for simple floating point parameters.  Note that the actual parameter values and 
    most attributes are actaully held by the model - the parameter object just holds an index into the
    model (why did I do this?)
    '''

    def __init__(self,model,code,name,index,factor=1.0,format='{0}'):
        '''
        Create a parameter definition

        Parameters hold two version of a value, one used in fitting (fitValue), and one 
        value presented to the user and stored in the XML file (value).  These may differ 
        by a scale factor (only the fitValue is actually stored).

        Arguments:
            model     The function model using the parameter
            code      The code for the parameter, used in XML
            name      The name of the parameter
            index     The index of the parameter in the fucmtion definition
            factor    A factor by which the parameter is multiplied for converting to XML
            format    A format string for the parameter

        '''

        self._model=model
        self._code=code
        self._name=name
        self._index=index
        self._factor=factor
        self._format=format
        self._error=0.0
        self._calcdate=None
        self._isLinear=False
        self._covarIndex=-1
        self._saved=None

    def code( self ):
        '''
        Return the code for the parameter
        '''
        return self._code

    def name( self ):
        '''
        Return the name of the parameter
        '''
        return self._name

    def fixed( self ):
        '''
        Return the fixed/calculate option for the parameter
        '''
        return self._model._fixed[self._index]

    def setFixed( self, fixed ):
        '''
        Select whether the parameter is to be fixed
        '''
        self._model._fixed[self._index]=bool(fixed)

    def setValue( self, valuestr ):
        '''
        Set the (user) value for the parameter.  Also sets the fit value, and
        clears the covariance index
        '''
        value=float(valuestr)
        value /= self._factor
        self.setFitValue(value)
        self._covarIndex=-1
        self._calcdate=None

    def getValue( self ):
        '''
        Get the user value
        '''
        value=self.fitValue()
        value *= self._factor
        return self._format.format(value)

    def getError( self ):
        '''
        Return the error (scaled to match the "user" value. 

        The error is from the last fit in which the parameter was calculated
        '''
        if self._error is None:
            return ''
        return self._format.format(self._error*self._factor)

    def covarIndex( self ):
        '''
        Return the index of the parameter in the covariance matrix from the last fit.
        '''
        return self._covarIndex

    def calcDate( self ):
        '''
        Return the date of the last model fit in which the parameter was calculated.
        '''
        return self._calcdate

    def fitValue( self ):
        '''
        The value of the parameter as used internally
        '''
        return self._model._param[self._index]

    def setFitValue( self, value, index=-1, error=None ):
        '''
        Updates the parameter, optionally setting the error and covariance index
        '''
        self._model._param[self._index] = value
        # Assume that if error is provided it has been calculated, so reset calc date
        self._error=error
        if index >= 0 and error != None:
            self._covarIndex=index
            self._calcdate=datetime.now()

    def saveValue( self ):
        '''
        Saves the current value to allow for restoring the value if a fit fails.
        '''
        self._saved=[self.fitValue(),self._error,self._calcdate,self._covarIndex]

    def restoreValue( self ):
        '''
        Restores the saved value
        '''
        if self._saved is not None:
            value,self._error,self._calcdate,self._covarIndex = self._saved
            self.setFitValue(value)

    def toXmlValue( self ):
        return self.getValue()

    def fromXmlValue( self, xmlstr ):
        self.setValue( xmlstr )

    def xmlElement( self ):
        element=ElementTree.Element('parameter')
        element.set('code',self.code())
        element.set('value',self.toXmlValue())
        element.set('fit','no' if self.fixed() else 'yes')
        if self._error:
            element.set('error',self.getError())
        if self._calcdate != None:
            element.set('calc_date',self._calcdate.strftime(datetimeformat))
        if self._covarIndex >= 0:
            element.set('covariance_index',str(self._covarIndex))
        return element

    def loadXmlElement( self, element ):
        assert(element.tag=='parameter')
        assert(element.get('code')==self.code())
        self.fromXmlValue(element.get('value',self.toXmlValue()))
        self.setFixed(element.get('fit','yes') != 'yes')
        self._error=None
        self._calcdate=None
        self._covarIndex=-1
        error=element.get('error','')
        calcdate=element.get('calc_date','')
        covarIndex=element.get('covariance_index','')
        if error:
            self._error = float(error)/self._factor
        if calcdate:
            try:
                self._calcdate=datetime.strptime(calcdate,datetimeformat)
            except:
                pass
        if covarIndex:
            try:
                self._covarIndex=int(covarIndex)
            except:
                pass

    def __str__( self ):
        return self.getValue()

class offset_parameter( parameter ):
    '''
    Offset parameter - uses millimetres as the user representation of the offset.
    '''
    def __init__(self,model,code,name,index):
        '''
        Create an offset parameter.

        Args:
            model    The model in which the offset will be used
            code     The code for the parameter (unique to the component)
            name     The name used to describe the parameter
            index    The index of the parameter in the component

        '''
        parameter.__init__( self,model,code,name,index,factor=1000,format="{0:.1f}")
        self._isLinear=True
        
class velocity_parameter( parameter ):
    '''
    Offset parameter - uses mm/year as the user representation of the offset, m/day internally.
    '''
    def __init__(self,model,code,name,index):
        '''
        Create an velocity parameter.

        Args:
            model    The model in which the offset will be used
            code     The code for the parameter (unique to the component)
            name     The name used to describe the parameter
            index    The index of the parameter in the component

        '''
        parameter.__init__( self,model,code,name,index,factor=1000*365.25,format="{0:.3f}")
        self._isLinear=True


class date_parameter( parameter ):
    '''
    Date parameter - uses a formatted date string as a user representation, and a day number internally.
    '''
    def __init__(self,model,code,name,index,format="{0:.1f}"):
        '''
        Create an date parameter.

        Args:
            model    The model in which the offset will be used
            code     The code for the parameter (unique to the component)
            name     The name used to describe the parameter
            index    The index of the parameter in the component

        '''
        parameter.__init__(self,model,code,name,index)

    def setValue(self, valuestr ):
        self.setFitValue(asday(valuestr))

    def getValue(self):
        return fromday(self.fitValue()).strftime("%d-%m-%Y %H:%M")

    def getError( self ):
        return '' if self._error is None else "{0:.2f}".format(self._error)

    def toXmlValue(self):
        return fromday(self.fitValue()).strftime(datetimeformat)


# Functions used to build model

class base_function( object ):
    '''
    Abstract base class for the functions (components) used to build a coordinate prediction model

    Each derived class implements a function calc for calculating the component at a given time, 
    and defines the parameters used by the component.

    The component is initiallized with with the model to which it will belong, a
    reference date, and the number of parameters it is defined by.  The reference date is a 
    parameter, which applies for most models (but not the offset, velocity, and cyclic components.
    It is used as a reference date for sorting the components, and more most models is a 
    fittable parameter of the model.  It is implemented as the last parameter of the model.

    It may be enabled or disabled - if it is disabled then it is not used in any calculations
    or fitting.  It may also be fitted or not fitted - if it is not fitted then it will not 
    generally be used in fitting the model.

    Each component parameter may also be fitted or not fitted.
    '''

    def __init__( self, model, date, nparam, params=None ):
        '''
        Create a base function - this is called by the specific component constructors.

        Args:
            model    The model to which the component will be attached
            date     The significant date of the component - for most component types this
                     is a parameter of the component (eg the date of an offset)
            nparam   The number of parameters in the model
            params   Optionally the parameters of the model

        '''
        if not date:
            date = model.refdate
        self.model = model
        self.parameters=[]
        self._event=''
        self._nparam = nparam+1
        self._fixed = [False] * (nparam+1)
        self._param = [0.0] * (nparam+1)
        self._param[nparam] = asday(date)
        self._fixed[nparam] = True
        self._funcFixed = True
        self._enabled = True
        if params:
            self._param[:nparam] = params[:nparam]

    def setModel( self, model ):
        '''
        Sets the parent model
        '''
        self.model = model

    def eventDate( self ):
        '''
        Returns the significant date of the component as a datetime object
        '''
        return fromday(self._param[-1])

    def _dateOffset( self, dates ):
        d = days_array(dates)-self._param[-1]
        return d.reshape(d.size,1)

    def setComponent( self, i, value, fixed=False ):
        '''
        Set the value of the i'th component and whether it is to be fixed or calculated.
        '''
        assert i >= 0 and i < self._nparam-1
        self._param[i] = float(value)
        self._fixed[i] = bool(fixed)

    def fixed(self):
        '''
        Returns the fixed status of the component
        '''
        return self._funcFixed

    def setFixed( self, fixed):
        '''
        Sets the component to be fixed (true) or calculated (false)
        '''
        self._funcFixed = fixed

    def enabled( self ):
        '''
        Returns the enabled status of the component (true or false)
        '''
        return self._enabled

    def setEnabled( self, enabled ):
        '''
        Sets the enabled status of the component (true or false)
        '''
        self._enabled = enabled

    def setEventName( self, name, detail='' ):
        '''
        Sets the name of the event the component represents
        '''
        self._event=name

    def eventName( self ):
        '''
        Returns the name of the event the component represents
        '''
        return self._event

    def eventDetail( self ):
        '''
        Returns more detailed information about the event the component represents
        '''
        return self._eventDetail

    def componentType( self ):
        '''
        Returns the type of the component.
        '''
        return type(self).__name__

    def xmlElement( self ):
        '''
        Function to returns an XML element representing the component.
        '''
        element=ElementTree.Element('component')
        element.set('type',type(self).__name__)
        if self._event:
            element.set('name',self._event)
        if not self._enabled:
            element.set('exclude','yes')
        element.set('fit','no' if self._funcFixed else 'yes')
        for p in self.parameters:
            element.append(p.xmlElement())
        return element

    @staticmethod
    def fromXmlElement( model, element ):
        '''
        Static function to reconstruct a componet from an XmlElement.
        '''
        classname = element.get('type','')
        if not classname:
            raise ValueError('Missing component type')
        def getsubclass( base, classname ):
            for c in base.__subclasses__():
                if c.__name__==classname:
                    return c
                subc=getsubclass(c,classname)
                if subc:
                    return subc
        cclass=getsubclass(base_function,classname)
        if not cclass:
            raise ValueError('Invalid station model type '+classname)
        component=cclass(model)
        component._event=element.get('name','')
        component._enabled=element.get('exclude','No').lower() != 'yes'
        component._funcFixed=element.get('fit','No').lower() != 'yes'
        for pelement in element:
            code=pelement.get('code','')
            for p in component.parameters:
                if p.code() == code:
                    p.loadXmlElement(pelement)
                    break
        return component

    def paramValues( self ):
        '''
        Return an array of parameter values for the model parameters.
        '''
        return [p.fitValue() for p in self.parameters]

    def paramStr( self ):
        '''
        Return a compact string representation of the component
        '''
        return (self.eventDate().strftime('%d-%m-%Y') + ' [' +
            ', '.join([str(p) for p in self._param[:-1]]) + ']')

    def __str__( self ):
        strval = self._event + ': ' if self._event else ''
        return strval + self.paramStr()
    
class offset( base_function ):
    '''
    A simple constant offset componet of the model (simplistically represents 
    the mean offset from the model reference coordinate)
    '''

    def __init__( self, model, offset=None ):
        base_function.__init__( self, model, refdate, 3, offset )
        self.parameters=[
            offset_parameter(self,'de_mm','East (mm)',0),
            offset_parameter(self,'dn_mm','North (mm)',1),
            offset_parameter(self,'du_mm','Up (mm)',2),
            ]
        self._funcFixed = False

    def calc( self, date ):
        return np.array([self._param[:3]]*len(date))

    def paramStr( self ):
        return type(self).__name__+ '  [{0:.1f}, {1:.1f}, {2:.1f}] mm'.format(
            *[x*1000 for x in self._param[:3]])

class velocity( base_function ):
    '''
    A simple constant velocity componet of the model (simplistically represents 
    the mean velocity of the mark)
    '''

    def __init__( self, model, velocity=None ):
        base_function.__init__( self, model, None, 3, velocity )
        self.parameters=[
            velocity_parameter(self,'ve_mmpy','East (mm/yr)',0),
            velocity_parameter(self,'vn_mmpy','North (mm/yr)',1),
            velocity_parameter(self,'vu_mmpy','Up (mm/yr)',2),
            ]
        self._funcFixed = False

    def calc( self, date ):
        y=self._dateOffset(date)
        return y.dot([self._param[:3]])

    def paramStr( self ):
        return type(self).__name__+ '  [{0:.2f}, {1:.2f}, {2:.2f}] mm/year'.format(
            *[x*daysperyear*1000 for x in self._param[:3]])

class velocity_change( base_function ):
    '''
    A velocity change component of the model.  This is added to the velocity
    for all dates after the component date.
    '''

    def __init__( self, model, date=refdate, velocity_change=None ):
        base_function.__init__( self, model, date, 3, velocity_change )
        self.parameters=[
            date_parameter(self,'date','Date of change',3),
            velocity_parameter(self,'ve_mmpy','East change (mm/yr)',0),
            velocity_parameter(self,'vn_mmpy','North change (mm/yr)',1),
            velocity_parameter(self,'vu_mmpy','Up change (mm/yr)',2),
            ]

    def calc( self, date ):
        y=np.maximum(self._dateOffset(date),0.0)
        return y.dot([self._param[:3]])

    def paramStr( self ):
        f=daysperyear*1000
        return type(self).__name__+ '  [{0:.2f}, {1:.2f}, {2:.2f}] mm/year at {3:%d-%m-%Y}'.format(
            self._param[0]*f, self._param[1]*f, self._param[2]*f, self.eventDate())

class cyclic( base_function ):
    '''
    Base class for the cyclic (annual, semi-annual) components
    '''

    def __init__( self, model, frequency, cosp=None, sinp=None ):
        if cosp and sinp:
            params = [cosp[0],cosp[1],cosp[2],sinp[0],sinp[1],sinp[2]]
        else:
            params=None
        base_function.__init__( self, model, refdate, 6, params )
        self._frequency=frequency
        self.parameters=[
            offset_parameter(self,'ecos_mm','East cosine (mm)',0),
            offset_parameter(self,'esin_mm','East sine (mm)',3),
            offset_parameter(self,'ncos_mm','North cosine (mm)',1),
            offset_parameter(self,'nsin_mm','North sine (mm)',4),
            offset_parameter(self,'ucos_mm','Up cosine (mm)',2),
            offset_parameter(self,'usin_mm','Up sine (mm)',5),
            ]

    def calc( self, date ):
        y=(self._dateOffset(date)/daysperyear)*math.pi*2*self._frequency;
        return np.cos(y).dot([self._param[:3]])+np.sin(y).dot([self._param[3:6]])

    def setComponent( self, i, issine, value, fixed=False ):
        '''
        Overwrite the setComponent function to allow setting the sine or cosine using 
        an index for the ordinate (E,N,U) and a boolean for sine (true) or cosine (false)
        component.
        '''
        assert type(issine) == bool, repr(issine)
        assert i >= 0 and i < 3
        if issine:
            i = i+3
        base_function.setComponent(self,i,value,fixed)

    def paramStr( self ):
        f=daysperyear*1000
        return (type(self).__name__+ 
                '  [{0:.2f}*cos+{3:.2f}*sin, {1:.2f}*cos+{4:.2f}*sin, {2:.2f}*cos+{5:.2f}*sin] mm'.format(
                    *[x*1000 for x in self._param[:6]])+
                ' frequency {0:.0f}/year'.format(self._frequency)
               )

class annual( cyclic ):
    '''
    Realisation of an annual cyclic component
    '''

    def __init__(self,model, cosp=None,sinp=None):
        cyclic.__init__(self,model, 1.0,cosp,sinp)

class semiannual( cyclic ):
    '''
    Realisation of an semi-annual cyclic component
    '''

    def __init__(self,model, cosp=None,sinp=None):
        cyclic.__init__(self,model, 2.0,cosp,sinp)

class equipment_offset( base_function ):
    '''
    A simple offset applying to all observations after a specified date.  Semantically
    represents a change due to equipment changes.
    '''

    def __init__( self, model, date=refdate, offset=None ):
        base_function.__init__( self, model, date, 3, offset )
        self.parameters=[
            date_parameter(self,'date','Date of offset',3),
            offset_parameter(self,'de_mm','East (mm)',0),
            offset_parameter(self,'dn_mm','North (mm)',1),
            offset_parameter(self,'du_mm','Up (mm)',2),
            ]

    def calc( self, date ):
        y=np.where(self._dateOffset(date)>0,1,0)
        return y.dot([self._param[:3]])

    def paramStr( self ):
        f=1000;
        return type(self).__name__+ '  [{0:.2f}, {1:.2f}, {2:.2f}] mm at {3:%d-%m-%Y}'.format(
            self._param[0]*f, self._param[1]*f, self._param[2]*f, self.eventDate())

class tectonic_offset( equipment_offset ):
    '''
    A simple offset applying to all observations after a specified date.  Semantically
    represents a change due to a tectonic event.
    '''

    def __init__( self, model, date=refdate, offset=None ):
        equipment_offset.__init__( self, model, date, offset )
        self.parameters[0]=date_parameter(self,'date','Date of event',3)

class slow_slip( base_function ):
    '''
    A slow slip event represented by an error function.
    '''

    def __init__( self, model, date=refdate, duration=10.0/3.92, offset=None ):
        base_function.__init__( self, model, date, 4 )
        self.setDuration(duration)
        if offset:
            self._param[:3]=offset[:3]
        self.parameters=[
            date_parameter(self,'mid_date','Central date of event',4),
            parameter(self,'duration_days','Duration of event (95% of slip) (days)',3,factor=3.92,format='{0:.1f}'),
            offset_parameter(self,'de_mm','East (mm)',0),
            offset_parameter(self,'dn_mm','North (mm)',1),
            offset_parameter(self,'du_mm','Up (mm)',2),
            ]

    def setDuration( self, duration, fixed=True ):
        self.setComponent(3,duration,fixed)

    def calc( self, date ):
        y=0.5*(1+erf(self._dateOffset(date)/self._param[3]))
        # y=np.minimum(1,np.maximum(0,self._dateOffset(date)/self._param[3]))
        return y.dot([self._param[:3]])
    
    def paramStr( self ):
        f=1000;
        start=self.eventDate()
        return type(self).__name__+ '  [{0:.2f}, {1:.2f}, {2:.2f}] mm centred {3:%d-%m-%Y} 95% width {4:.2f} days'.format(
            self._param[0]*f, self._param[1]*f, self._param[2]*f, start, 3.92*self._param[3])

class slow_slip_ramp( base_function ):
    '''
    A slow slip event represented by a simple ramp function between the start and end date
    '''

    def __init__( self, model, date=refdate, end_date=None, offset=None ):
        base_function.__init__( self, model, date, 4 )
        if end_date is None:
            end_date=self._param[4]+10.0
        self.setEndDate(end_date)
        if offset:
            self._param[:3]=offset[:3]
        self.parameters=[
            date_parameter(self,'start_date','Start date of event',4),
            date_parameter(self,'end_date','End date of event',3),
            offset_parameter(self,'de_mm','East (mm)',0),
            offset_parameter(self,'dn_mm','North (mm)',1),
            offset_parameter(self,'du_mm','Up (mm)',2),
            ]

    def setEndDate( self, end_date, fixed=True ):
        self.setComponent(3,asday(end_date),fixed)

    def calc( self, date ):
        duration=self._param[3]-self._param[4]
        y=self._dateOffset(date)/duration
        y=np.minimum(1.0,np.maximum(y,0.0))
        return y.dot([self._param[:3]])

    def paramStr( self ):
        f=1000;
        start=self.eventDate()
        return type(self).__name__+ '  [{0:.2f}, {1:.2f}, {2:.2f}] mm between {3:%d-%m-%Y} and {4:%d-%m-%Y}'.format(
            self._param[0]*f, self._param[1]*f, self._param[2]*f, start, fromday(self._param[3]))

class exponential_decay( base_function ):
    '''
    A component representing an exponential decaying rate of change from a given start date
    '''

    def __init__( self, model, date=refdate, decay=10.0/math.log(2.0), offset=None ):
        base_function.__init__( self, model, date, 4 )
        self.setDuration(decay)
        if offset:
            self._param[:3]=offset[:3]
        self.parameters=[
            date_parameter(self,'date','Start date of event',4),
            parameter(self,'halflife_days','Half life of event (days)',3,factor=math.log(2.0),format='{0:.1f}'),
            offset_parameter(self,'de_mm','East (mm)',0),
            offset_parameter(self,'dn_mm','North (mm)',1),
            offset_parameter(self,'du_mm','Up (mm)',2),
            ]

    def calc( self, date ):
        y=self._dateOffset(date)/self._param[3]
        y=np.maximum(0,1.0-np.exp(-y))
        return y.dot([self._param[:3]])

    def setDuration( self, decay, fixed=True ):
        self.setComponent(3,decay,fixed)

    def paramStr( self ):
        f=1000;
        start=self.eventDate()
        return type(self).__name__+ '  [{0:.2f}, {1:.2f}, {2:.2f}] mm exponential decay from {3:%d-%m-%Y} half life {4:.1f} days'.format(
            self._param[0]*f, self._param[1]*f, self._param[2]*f, start, self._param[3]*math.log(2.0))

class exclude_obs( object ):
    '''
    Class representing an observation excluded at a specific date
    '''

    def __init__( self, date, comment="", index=-1 ):
        self.date=date
        self.day=asday(date)
        self.index=index
        self.comment=comment

    def xmlElement(self):
        element=ElementTree.Element('exclude')
        element.set('date',self.date.strftime('%Y-%m-%dT%H:%M:%S'))
        if self.comment:
            element.set('comment',self.comment)
        return element

    @staticmethod
    def fromXmlElement( element ):
        assert element.tag == 'exclude'
        date = datetime.strptime(element.get('date'),'%Y-%m-%dT%H:%M:%S')
        comment = element.get('comment','')
        return exclude_obs( date, comment )

    
class model( object ):
    '''
    Class representing a station prediction model.  

    The class contains a coordinate prediction model defining the predicted offset
    of the mark from its reference coordinate at any given date.  The model has a number of 
    components, as a minimum it always has an offset, velocity, annual, and semiannual term.
    Additional components can be added and removed with the addComponent and removeComponent
    functions.

    The model can be evaluated at a set of dates using the calc function, returning either 
    XYZ coordinates, or ENU offsets from the reference point.

    Observations can be loaded into the model from a time series file using the 
    loadTimeSeries function. The model components can then be fitted to the time series using 
    the fit and fitAllLinear functions.  The loaded time series can be extracted with the 
    getObs function.

    The class can also hold a list of dates for which the time series data are rejected from
    fitting calculations, and a set of outages - periods for which the time series contains
    no data.

    The class can be persisted to an XML file.  The time series data is not persisted - that 
    must be reloaded each time the model is instantiated if it is required.
    '''

    BasicComponents=[offset,velocity,annual,semiannual]

    def __init__(self, station=None, xyz=None, filename=None, loadfile=True ):
        self.refdate = refdate
        self.versiondate = datetime.now()
        self.setStation(station,xyz)
        self.filename=None
        self.dates=None
        self.obs=None
        self.useobs=None
        self.excluded=[]
        self.components=[c(self) for c in self.BasicComponents]
        self.covariance=None
        if filename and loadfile:
            self.load(filename)
        self.filename=filename
        self.saved=self.toXmlString()

    def setStation( self, station, xyz, site='', priority=1 ):
        '''
        Reset the station code for the model. 

        Args:
            station    The station code
            xyz        The station location (reference coordinate)
            site       The station site (code for a group of stations)
            priority   The priority of the station within the site (lowest priority
                       is preferred station). Should be an integer
        '''
        self.station = station or ''
        self.xyz=xyz
        self.site=site or station
        self.priority=int(priority)
        grs80 = ellipsoid.grs80
        if xyz is not None:
            lon,lat,h=grs80.geodetic(xyz)
            self.enu_axes=grs80.enu_axes(lon,lat)

    def copy( self, filename=None ):
        '''
        Return a copy of the current model.  
        
        Optionally sets a filename that will be used for saving the model.
        The copy includes the model components and the excluded observations, 
        but not other information (eg loaded time series, outages).

        Args:
            filename    An optional filename for the model

        Returns:
            A new stn_pred_model.model object with the same parameters set
            as the current model
        '''
        xml=self.toXmlElement()
        copy=model(station=self.station,filename=filename,loadfile=False)
        copy.loadFromXml(xml,'copy')
        return copy

    def toXmlString( self ):
        return ElementTree.tostring(self.toXmlElement())

    def updateXmlAttributes( self, root ):
        root.set('code',self.station)
        root.set('site',self.site)
        root.set('priority',str(self.priority))
        root.set('version_date',self.versiondate.strftime(datetimeformat))

    def toXmlElement( self ):
        root=ElementTree.Element(stn_tag)
        self.updateXmlAttributes( root )
        
        spm=ElementTree.Element(cpm_tag)
        spm.set('ref_date',self.refdate.strftime(datetimeformat))
        if self.xyz != None:
            spm.set('X0',str(self.xyz[0]))
            spm.set('Y0',str(self.xyz[1]))
            spm.set('Z0',str(self.xyz[2]))
        components=ElementTree.Element('components')
        for c in self.components:
            components.append(c.xmlElement())
        spm.append(components)
        if self.covariance is not None:
            c=self.covariance
            size=c.shape[0]
            covar=ElementTree.Element('covariance')
            covar.set('size',str(size))
            for i in range(size):
                for j in range(i+1):
                    el=ElementTree.Element('element')
                    el.set('row',str(i))
                    el.set('col',str(j))
                    el.set('value',str(c[i,j]))
                    covar.append(el)
            spm.append(covar)
        if self.excluded:
            excluded=ElementTree.Element('excluded')
            for e in self.excluded:
                excluded.append(e.xmlElement())
            spm.append(excluded)
        root.append(spm)
        return root

    def changed( self ):
        xmlstr=self.toXmlString()
        return xmlstr != self.saved

    def getFilename( self, filename=None ):
        if not filename:
            filename=self.filename
        if filename and self.station:
            filename=filename.replace('{code}',self.station)
        return filename

    def readStationXmlFile(self,filename):
        if not os.path.exists(filename):
            raise RuntimeError('Station file '+filename+' does not exist')
        with open(filename) as mf:
            xmlstr=mf.read()
            xmlstr=re.sub(r'\>\s+\<','><',xmlstr)
        root=ElementTree.fromstring(xmlstr)
        if root.tag != stn_tag:
            raise ValueError(filename+' is not a station file')
        code = root.get('code','')
        if not code:
            raise ValueError(filename+' does not specify a station code')
        if self.station and code != self.station:
            raise ValueError(filename+' is not for station '+self.station)
        return root

    def save( self, filename=None, updateAvailability=False ):
        filename = self.getFilename( filename )
        if not filename:
            raise ValueError('No file name specified for saving station prediction model')

        self.versiondate=datetime.now()
        root = self.toXmlElement()
        if os.path.exists(filename):
            oldroot=self.readStationXmlFile(filename)
            oldspm=oldroot.find(cpm_tag)
            if oldspm is not None:
                 oldroot.remove(oldspm)
            newspm=root.find(cpm_tag)
            oldroot.append(newspm)
            root=oldroot

        self.updateXmlAttributes( root )

        if updateAvailability:
            self.updateAvailability(root)

        with open(filename,'w') as f:
            xmlstr=ElementTree.tostring(root)
            pxmlstr=minidom.parseString(xmlstr).toprettyxml(indent='  ')
            f.write(pxmlstr)
            self.saved = self.toXmlString()

    def load( self, filename ):
        filename = self.getFilename( filename )
        if not filename:
            raise ValueError('No file name specified for loading station prediction model')
        if not os.path.exists(filename):
            raise RuntimeError('Station prediction model file '+filename+' does not exist')
        
        root=self.readStationXmlFile(filename)
        self.loadFromXml(root,filename)
        self.filename=filename
        self.saved=self.toXmlString()

    def loadFromXmlString( self, xmlstr, source="XML string" ):
        root=ElementTree.fromstring(xmlstr)
        self.loadFromXml(root)

    def loadFromXml( self, root, source="XML string" ):
        #tree=ElementTree.ElementTree(ElementTree.fromstring(xmlstr))
        #root=tree.getroot()
        if root.tag != stn_tag:
            raise ValueError(source+' is not station station')
        code = root.get('code','')
        if not code:
            raise ValueError(source+' does not specify a station code')
        site = root.get('site',code)
        priority=int(root.get('priority','1'))
        if self.station and code != self.station:
            raise ValueError(source+' is not for station '+self.station)

        vrdt=root.get('version_date')
        if vrdt:
            self.versiondate=datetime.strptime(vrdt,datetimeformat)
            
        spm=root.find(cpm_tag)
        rfdt=spm.get('ref_date')
        if rfdt:
            self.refdate=datetime.strptime(rfdt,datetimeformat)

        xyz=[0,0,0]
        for i,axis in enumerate(('X0','Y0','Z0')):
            try:
                xyz[i]=float(spm.get(axis,''))
            except:
                raise ValueError(source+' does not define an '+axis+' value')
        self.setStation(code,xyz)

        components=[]
        comproot=spm.find('components')
        for c in comproot:
            components.append(base_function.fromXmlElement(self,c))
        self.components=components
        self.addBasicComponents()

        self.covariance=None
        covar=spm.find('covariance')
        if covar is not None:
            size=int(covar.get('size'))
            self.covariance=np.zeros((size,size))
            for el in covar:
                if el.tag == 'element':
                    r=int(el.get('row'))
                    c=int(el.get('col'))
                    v=float(el.get('value'))
                    self.covariance[r,c]=v
                    self.covariance[c,r]=v

        self.excluded=[]
        excluded=spm.find('excluded')
        if excluded is not None:
            for e in excluded:
                self.excluded.append(exclude_obs.fromXmlElement(e))
        self.setExcludedObs()

    def sortComponents( self ):
        '''
        Sorts components.  The basic components are sorted to the front
        of the list, then the others by their reference date (assumed to
        be the last parameter)
        '''
        bc=self.BasicComponents
        def keyf(comp):
            ct=type(comp)
            return (
                bc.index(ct) if ct in bc else len(bc),
                int(comp._param[-1]),
                ct.__name__,
                comp._param[-1]
                )
        self.components.sort(key=keyf)

    def addBasicComponents( self ):
        '''
        Adds the basic components all models use (ie offset, velocity, annual, and 
        semi-annual components)
        '''
        ctypes = [type(c) for c in self.components]
        for btype in self.BasicComponents:
            if btype not in ctypes:
                component=btype(self)
                if btype in (annual,semiannual):
                    component.setEnabled(False)
                self.components.append(component)
        self.sortComponents()

    def addComponent( self, component ):
        '''
        Adds a new component to the model

        The component must be initiallized to reference this model
        '''
        if component in self.components:
            return
        if component.model != self:
            raise RuntimeError('Cannot add component for a different model')
        self.components.append(component)
        self.sortComponents()

    def removeComponent( self, component ):
        '''
        Removes a component from the model
        '''
        if component in self.components:
            self.components.remove(component)
            component.model=None

    def calc( self, dates, enu=True ):
        '''
        Calculate the model at one or more dates. 

        Args:
            dates   Either a single date or a list of dates
            enu     If true then returns the E,N,U components relative to the 
                    reference coordinate.  Otherwise returns geocentric
                    X,Y,Z values

        Returns:
            Either a single coordinate or an array of coordinates
        '''
        single=not isinstance(dates,list) and not isinstance(dates,np.ndarray)
        dates=days_array(dates)
        value=np.zeros((len(dates),3))
        for m in self.components:
            if m.enabled():
                value += m.calc(dates)
        if not enu:
            value=self.xyz+value.dot(self.enu_axes)

        return value[0] if single else value

    def setUseObs( self, index, comment=None, use=True, toggle=False ):
        '''
        Sets an observation to be used or not used for model fitting.

        If the observation is set to not used, then its date is added to
        the list of excluded dates for the model.

        Args:
            index   The index of the observation in the time series.
            comment A comment describing why the observation is not used
            use     True or false to use or not use the observation
            toggle  If true then the current value is toggled
    
        '''
        if toggle:
            use = not self.useobs[index]
        elif use == self.useobs[index]:
            return
        self.useobs[index]=use
        found = False
        if use:
            for e in self.excluded:
                if e.index==index:
                    self.excluded.remove(e)
                    break
        else:
            date=self.dates[index]
            self.excluded.append(exclude_obs(date,comment,index=index))
            self.excluded.sort( key=lambda x: x.day )

    def setExcludedObs( self ):
        '''
        Used to initially link the excluded obs list to the time series data
        '''
        if self.dates is not None:
            self.useobs=np.array([True]*len(self.dates))
            days=[asday(d) for d in self.dates]
            for e in self.excluded:
                e.index=-1
                for i,d in enumerate(days):
                    if abs(d-e.day) < 0.5:
                        e.index=i
                        self.useobs[i]=False
                        break

    def loadTimeSeries( self, filename, transform=None ):
        '''
        Loads a time series to be analysed with the model

        Assumes the time series file has columns name, epoch, x, y, z and is 
        whitespace delimited.  The filename can include {code}, which will be replaced
        with the station code associated with the model.

        Flags observations as excluded if they match the dates stored with the model.
        '''

        if self.station and '{code}' in filename:
            filename = filename.replace('{code}',self.station)
        with open(filename) as f:
            fields=f.readline().split()

            if fields[:5] != ('name epoch x y z'.split()):
                raise RuntimeError('Station file '+filename+' doesn\'t have the correct fields')

            obs=[]
            tsobs=namedtuple('tsobs','epoch enu')
            for l in f:
                parts=l.split()
                if len(parts) < 5:
                    continue
                if self.station==None:
                    self.station=parts[0].upper()
                if parts[0].upper() != self.station:
                    continue
                epoch = datetime.strptime(parts[1],datetimeformat)
                xyz = np.array([float(p) for p in parts[2:5]])
                if transform:
                    xyz=transform(xyz,epoch)
                if self.xyz==None:
                    self.setStation(self.station,xyz)
                enu=self.enu_axes.dot(xyz-self.xyz)
                obs.append(tsobs(epoch,enu))

            obs.sort(key=lambda o: o.epoch)
            self.dates=np.array([o.epoch for o in obs])
            self.enu=np.array([o.enu for o in obs])
            self.setExcludedObs()

    def getObs( self ):
        '''
        Returns the currently loaded time series.

        Returns three values:
            dates    An array of dates for the time series
            enu      An array of east,north,up values relative to the
                     model reference coordinate
            useobs   A boolean array flagging the usage of the observation

        '''
        return self.dates,self.enu,self.useobs

    @staticmethod
    def robustStandardError(obs):
        '''
        Estimate the standard error for components of a time series.

        Standard errors are estimated using the differences between consecutive elements
        of the time series.  The 95 percentile of the absoluted differences is used as a
        robust estimator for the 95% cumulative probability (two tailed), so is scaled by 
        1.96 to get standard error.
        '''
        errors=[0]*3
        for axis in range(3):
            diffs=np.abs(obs[1:,axis]-obs[:-1,axis])
            # Note sqrt(2) accounts for the fact that these are differences
            # between two observations
            se=np.percentile(diffs,95.0)/(1.96*np.sqrt(2))
            errors[axis]=se
        return errors

    def clearCovariance( self ):
        '''
        Empties the covariance matrix - this is recomputed each time a fit is calculated
        '''
        self.covariance=None
        for m in self.components:
            for p in m.parameters:
                p._covarIndex=-1

    def fit( self ):
        '''
        Fit all flagged parameters of all flagged and enabled models
        (ie with fit flag set)
        '''
        # Form a list of parameters that we are fitting
        fit_params=[]
        for m in self.components:
            if not m.enabled():
                continue
            if not m.fixed():
                for p in m.parameters:
                    if not p.fixed():
                        fit_params.append(p)
        return self.fitParams( fit_params )

    def fitAllLinear( self, compo = None ):
        '''
        Fit all linear parameters of all enabled models. Allow non linear fitting for one component at a time
        '''
        # Form a list of parameters that we are fitting
        fit_params=[]
        for m in self.components:
            if not m.enabled():
                continue
            for p in m.parameters:
                if (p._isLinear or m == compo) and not p.fixed() :
                    fit_params.append(p)
        return self.fitParams( fit_params )

    def fitParams( self, fit_params ):
        '''
        Fit a selected set of parameters.  

        Args:
            fit_params   An array of parameters (of model functions)
                         that are to be fitted

        Returns:
            ok           Boolean success/failure status
            mesg         Status message

        '''
        if not fit_params:
            return True, 'Nothing to fit'

        # Build a set of fitted parameters, and for each 
        # parameter save it's currentvalue
        fitting=set()
        for p in fit_params:
            p.saveValue()
            fitting.add(p._model)

        # Start values the initial values for the (potentially) non-linear fit
        start_values = [p.fitValue() for p in fit_params]
        # Determine standard errors based on differences between obs
        # Used to weight observations in fit
        se = np.array([self.robustStandardError(self.enu)])

        # Form the working arrays of observation dates, ENU coordiantes,
        # and usage flag
        dates=days_array(self.dates)
        res=self.enu
        useobs=self.useobs
    
        # Correct the obs for the components we are not fitting
        first=True
        for m in self.components:
            if not m.enabled():
                continue
            if m in fitting:
                continue
            if first:
                res=self.enu.copy()
                first=False
            res -= m.calc(dates)

        # Function to calc the residuals.  This is used by the fitting routine.  It
        # returns a one dimensional array of residuals containing each component for 
        # each observation.  The components are scaled by the estimated standard errors
        # for the component (ie individually for E, N, and U)

        def calcres(params):
            # Assign param values
            for p,v in zip(fit_params,params):
                p.setFitValue(v)
            # Calculate the residual vector
            vres=res.copy()
            for m in fitting:
                vres -= m.calc(dates)
            vres /= se
            vres[~useobs]=[0,0,0]
            return vres.reshape((vres.size,))

        # Use leastsq function to do the fit ...
        x, covx, info, mesg, ier = leastsq(calcres,start_values,full_output=1)

        ok = True
        if ier in [1,2,3,4] and covx is not None:
            mesg = 'Model fitted successfully'
            self.clearCovariance()
            for i,p in enumerate(fit_params):
                v=x[i]
                error=np.sqrt(covx[i,i])
                p.setFitValue(v,i,error)
                self.covariance=covx
        else:
            # It not successful, then restore the original values
            mesg = 'Model not fitted: '+str(mesg)
            ok = False
            for p in fit_params:
                p.restoreValue()

        return ok, mesg

    def updateAvailability(self,root):
        '''
        Update the xml object with the outages in the time series

        The outages are potentially used to determine when data from the station
        might (or might not) be available.  These are stored in the model XML file.
        '''
        if self.dates is None or len(self.dates) < 1:
            return
        firstobs=self.dates[0].date()
        outages=ElementTree.Element(outages_tag)
        oneday=timedelta(days=1)
        nextdate=firstobs+oneday
        startofday=timetype(0,0,0)
        endofday=timetype(23,59,59)
        for epoch in self.dates[1:]:
            obsdate=epoch.date()
            if obsdate > nextdate:
                start=datetime.combine(nextdate,startofday)
                end=datetime.combine(obsdate-oneday,endofday)
                outage=ElementTree.Element(outage_tag)
                outage.set('start',start.strftime(datetimeformat))
                outage.set('end',end.strftime(datetimeformat))
                outages.append(outage)
            nextdate=obsdate+oneday

        root.set('start_date',datetime.combine(firstobs,startofday).strftime(datetimeformat))
        oldoutages=root.find(outages_tag)
        if oldoutages is not None:
            root.remove(oldoutages)
        root.append(outages)

    def __str__( self ):
        '''
        Print a readable description of the model.
        '''
        descr=['Station: '+str(self.station)]
        descr.extend([str(m) for m in self.components if m.enabled()])
        # descr.extend([str(self.events[k]) for k in sorted(self.events.keys())])
        return '\n    '.join(descr)


    def readGnsFiles( self, filename ):
        '''
        Loads a GNS model file, reading three components, E,N, and U

        Expects a file name with placeholder {enu} which will be substituted with e, n, and u
        to find the  3 files required.  eg PYGR_{enu}.out.

        This code was used to import the original station models generated by John Beavan.
        '''

        if '{code}' in filename and self.station:
            filename=filename.replace('{code}',self.station)

        if '{enu}' not in filename:
            raise ValueError('GNS stn prediction model filename must include {enu} placeholder: '+filename)

        events={}
        def _getEvent( type, date, *params ):
            key=type.__name__+str(int(date))+'_'.join("{:.1f}".format(x) for x in params)
            if key not in events:
                model=type(self,date,*params)
                events[key] = model
                self.components.append(model)
            return events[key]

        self.components=[c(self) for c in self.BasicComponents]

        axes=['e','n','u']
        mm = lambda x: float(x)/1000.0
        years = lambda x: float(x)*365.25
        tfixed = lambda x: not bool(int(x))
        def parseline( f, *types ):
            parts = f.readline().split()
            if len(parts) < len(types):
                raise ValueError
            return [t(p) for t,p in zip(types,parts)]
        
        for i,c in enumerate(axes):
            cfile = filename.replace('{enu}',c)
            with open(cfile) as f:
                header=f.readline()
                start_time=f.readline()
                end_time=f.readline()
                offset,offsetfixed=parseline(f,mm,tfixed)
                self.components[1].setComponent(i, *parseline(f,mm,tfixed))
                self.components[2].setComponent(i, False, *parseline(f,mm,tfixed))
                self.components[2].setComponent(i, True, *parseline(f,mm,tfixed))
                self.components[3].setComponent(i, False, *parseline(f,mm,tfixed))
                self.components[3].setComponent(i, True, *parseline(f,mm,tfixed))

                # Velocity change
                for nc in  range(*parseline(f,int)):
                    date,change,fixed=parseline(f,float,mm,tfixed)
                    _getEvent(velocity_change,date).setComponent(i,change,fixed)

                # Equipment offset
                for nc in  range(*parseline(f,int)):
                    date,change,fixed=parseline(f,float,mm,tfixed)
                    _getEvent(equipment_offset,date).setComponent(i,change,fixed)

                # Tectonic offset
                for nc in  range(*parseline(f,int)):
                    date,change,fixed=parseline(f,float,mm,tfixed)
                    _getEvent(tectonic_offset,date).setComponent(i,change,fixed)

                # Exponential
                for nc in  range(*parseline(f,int)):
                    date,duration,change,fixedd,fixedc=parseline(f,float,float,mm,tfixed,tfixed)
                    duration = 1.0/duration
                    component=_getEvent(exponential_decay,date,duration)
                    component.setDuration(duration,fixedd)
                    component.setComponent(i,change,fixedc)

                # Slow slip
                for nc in  range(*parseline(f,int)):
                    date,duration,change,fixedd,fixedc=parseline(f,float,float,mm,tfixed,tfixed)
                    duration = 1.0/duration
                    if duration < 0:
                        duration=-duration
                        change=-change
                    component=_getEvent(slow_slip,date,duration)
                    component.setDuration(duration,fixedd)
                    component.setComponent(i,change,fixedc)
                    # Slow slip calculated differently for GNS version - negative before 
                    # start of slip rather than 0...
                    offset -= change/2.0
            self.components[0].setComponent(i, offset, offsetfixed )
            # Decide whether we are using annual and semi-annual
            for ic in (2,3):
                c=self.components[ic]
                c.setEnabled(False)
                for p in c.parameters:
                    if not p.fixed() or p.fitValue() != 0.0:
                        c.setEnabled(True)
                        break
        self.sortComponents()

if __name__ == '__main__':
    '''
    This module can be used as a stand alone application to calculate the model
    at a set of dates.
    '''
    import argparse

    parser=argparse.ArgumentParser('Calculate a time series from a station prediction model')
    parser.add_argument('code',help='Code of station to calculate or the name of a model file')
    parser.add_argument('start_date',nargs='?',help='Start date for calculating (YYYY-MM-DD) or filename')
    parser.add_argument('end_date',nargs='?',help='End date for calculating (YYYY-MM-DD)')
    parser.add_argument('-f','--model-file',default='stations/{code}.xml',help="File name for model")
    parser.add_argument('-t','--timeseries-file',default='timeseries/{code}_igs08_xyz.dat',help="Time series file")
    parser.add_argument('-m','--model-dir',default='stations',help='Base directory for models (default stations)')
    parser.add_argument('-x','--calc-xyz',action='store_true',help='Calculate XYZ instead of enu')
    parser.add_argument('-i','--increment_days',type=int,help='Increment in days for calculation')
    parser.add_argument('-d','--debug-calcs',action='store_true',help='Print individual components')
    parser.add_argument('-u','--update_model',action='store_true',help='Load timeseries, fit linear components, update availability')

    args=parser.parse_args()

    model_file=args.model_dir+'/'+args.model_file
    code=args.code

    if os.path.isfile(code):
        model_file=code
        code=None

    try:
        spm=model(station=code,filename=model_file)

        if args.update_model:
            spm.loadTimeSeries(args.timeseries_file)
            spm.fitAllLinear()
            spm.save(updateAvailability=True)

        if args.start_date is None:
            print spm
        else:
            days=[]
            if args.end_date is None and os.path.exists(args.start_date):
                with open(args.start_date) as df:
                    for l in df:
                        m=re.match(r'^\s*(\d\d\d\d\-\d\d\-\d\d)\s*$',l)
                        if m:
                            days.append(datetime.strptime(m.group(1),'%Y-%m-%d'))
            else:
                sdate=datetime.strptime(args.start_date,'%Y-%m-%d')
                edate=datetime.strptime(args.end_date,'%Y-%m-%d') if args.end_date is not None else sdate
                incdays=args.increment_days
                if incdays <= 0:
                    raise RuntimeError('Date increment must be positive')
                tdel=timedelta(days=incdays)
                while sdate <= edate:
                    days.append(sdate)
                    sdate += tdel
            calcenu=not args.calc_xyz
            format="{0:%Y-%m-%d}\t{1:.#f}\t{2:.#f}\t{3:.#f}".replace('#','1' if calcenu else '4')
            for sdate in days:
                xyz=spm.calc(sdate,calcenu)
                if calcenu:
                    xyz *= 1000.0;
                print format.format(sdate,xyz[0],xyz[1],xyz[2])
                if args.debug_calcs:
                    for c in m.components:
                        if c.enabled():
                            xyz = c.calc(days_array(sdate))[0]*1000.0
                            print "{0}\t{1:.1f}\t{2:.1f}\t{3:.1f}".format(
                                c.componentType(),xyz[0],xyz[1],xyz[2])

    except:
        msg=str(sys.exc_info()[1])
        print "Error: "+msg
