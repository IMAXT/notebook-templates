from pathlib import Path
import holoviews as hv
from holoviews import dim, opts
from holoviews.operation.datashader import regrid, dynspread, datashade
from imaxt_image.image import TiffImage, zscale
from astropy.table import Table
import warnings
import numpy as np
import param
import panel as pn
from functools import partial
import panel.widgets as pnw

warnings.filterwarnings("ignore")

class zscale_filter(hv.Operation):

    normalize = param.Boolean(default=False)

    def _process(self, element, key=None):
        xs = element.dimension_values(0, expanded=False)
        ys = element.dimension_values(1, expanded=False)

        # setting flat=False will preserve the matrix shape
        data = element.dimension_values(2, flat=False)

        if self.p.normalize:
            dr = data.ravel()
            data = (data - dr.mean()) / dr.std() * 2 ** 16

        vmin, vmax = zscale(data.ravel())

        new_data = data.clip(vmin, vmax)

        label = element.label
        # make an exact copy of the element with all settings, just with different data and label:
        element = element.clone((xs, ys, new_data), label=label)
        return element
    
channel = pnw.IntSlider(name='Channel', value=1, start=1, end=15)
zslice = pnw.Select(name='Slice', options=[])

class CubeViewer(param.Parameterized):
    
    channel = param.Integer()
    zslice = param.ObjectSelector()
    
    def __init__(self, data_dir, imgtype='cubes', xaxis=None, yaxis=None, colorbar=True, toolbar='below', height=400):
        self.data_dir = Path(data_dir)
        self.imgtype = imgtype
        self.xaxis = xaxis
        self.yaxis = yaxis
        self.colorbar = colorbar
        self.toolbar = toolbar
        self.height = height
        self.slice_files = [f.name for f in (self.data_dir / self.imgtype).glob('*.tif')]
        
        self.param['channel'].bounds=(0, self.channels-1)
        self.param['zslice'].objects=self.slice_files
        self.param['zslice'].default=self.slice_files[0]
        super().__init__()
                
    @property
    def channels(self):
        image = TiffImage(self.data_dir / self.imgtype / self.slice_files[0])
        return image.shape[0]
    
    @param.depends('channel', 'zslice')
    def view(self):
        with TiffImage(self.data_dir / self.imgtype / self.zslice) as im:
            image = im.asarray()
            img = image[self.channel]
            xsize, ysize = img.shape
            img = hv.Image(img, bounds=(0, 0, ysize, xsize), vdims='Intensity')
            img = zscale_filter(img)
            return regrid(img).opts(plot={'Image': dict(frame_height=self.height,
                                                        colorbar=self.colorbar,
                                                         toolbar=self.toolbar,
                                                         xaxis=self.xaxis,
                                                         yaxis=self.yaxis,
                                                         aspect='equal')})
        
    def display(self, height=400):
        res = pn.Row(self.param, self.view)
        res.data_dir = self.data_dir
        return res
    
        
class Catalogue:
    def __init__(self, data_dir, filename, height=400):
        self.data_dir = Path(data_dir)
        self.height = height
        self.df = Table.read(self.data_dir / 'catalogue' / filename).to_pandas()
        
    def display(self):
        points = hv.Points((self.df['X'], self.df['Y'], np.log10(self.df['flux_25'])), vdims='z')
        res = points.opts(size=0.25) # dynspread(datashade(points))
        res.opts(aspect='equal').opts(opts.Points(color='z')).opts({'Points': {'frame_height': self.height}})
        return res

    @classmethod
    def from_view(cls, view):
        r = view.get_root()
        s = r.children[0].children[-1]
        f = view.get_root().children[-1].children[-1]
        filename = s.value.replace('.tif','.fits')
        return cls(view.data_dir, filename, height=f.frame_height)
        
        
    