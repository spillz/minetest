#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# This program is free software. It comes without any warranty, to
# the extent permitted by applicable law. You can redistribute it
# and/or modify it under the terms of the Do What The Fuck You Want
# To Public License, Version 2, as published by Sam Hocevar. See
# COPYING for more details.

# Made by Jogge, modified by celeron55
# 2011-05-29: j0gge: initial release
# 2011-05-30: celeron55: simultaneous support for sectors/sectors2, removed
# 2011-06-02: j0gge: command line parameters, coordinates, players, ...
# 2011-06-04: celeron55: added #!/usr/bin/python2 and converted \r\n to \n
#                        to make it easily executable on Linux
# 2011-07-30: WF: Support for content types extension, refactoring
# 2011-07-30: erlehmann: PEP 8 compliance.

# Requires Python Imaging Library: http://www.pythonware.com/products/pil/

# Some speed-up: ...lol, actually it slows it down.
#import psyco ; psyco.full()
#from psyco.classes import *

import zlib
import os
import string
import time
import getopt
import sys
import cStringIO
import traceback
import numpy
import itertools
from PIL import Image, ImageDraw, ImageFont, ImageColor


TRANSLATION_TABLE = {
    1: 0x800,  # CONTENT_GRASS
    4: 0x801,  # CONTENT_TREE
    5: 0x802,  # CONTENT_LEAVES
    6: 0x803,  # CONTENT_GRASS_FOOTSTEPS
    7: 0x804,  # CONTENT_MESE
    8: 0x805,  # CONTENT_MUD
    10: 0x806,  # CONTENT_CLOUD
    11: 0x807,  # CONTENT_COALSTONE
    12: 0x808,  # CONTENT_WOOD
    13: 0x809,  # CONTENT_SAND
    18: 0x80a,  # CONTENT_COBBLE
    19: 0x80b,  # CONTENT_STEEL
    20: 0x80c,  # CONTENT_GLASS
    22: 0x80d,  # CONTENT_MOSSYCOBBLE
    23: 0x80e,  # CONTENT_GRAVEL
    24: 0x80f,  # CONTENT_SANDSTONE
    25: 0x810,  # CONTENT_CACTUS
    26: 0x811,  # CONTENT_BRICK
    27: 0x812,  # CONTENT_CLAY
    28: 0x813,  # CONTENT_PAPYRUS
    29: 0x814}  # CONTENT_BOOKSHELF


def hex_to_int(h):
    i = int(h, 16)
    if(i > 2047):
        i -= 4096
    return i


def hex4_to_int(h):
    i = int(h, 16)
    if(i > 32767):
        i -= 65536
    return i


def int_to_hex3(i):
    if(i < 0):
        return "%03X" % (i + 4096)
    else:
        return "%03X" % i


def int_to_hex4(i):
    if(i < 0):
        return "%04X" % (i + 65536)
    else:
        return "%04X" % i


def getBlockAsInteger(p):
    return p[2]*16777216 + p[1]*4096 + p[0]

def unsignedToSigned(i, max_positive):
    if i < max_positive:
        return i
    else:
        return i - 2*max_positive

def getIntegerAsBlock(i):
    x = unsignedToSigned(i % 4096, 2048)
    i = int((i - x) / 4096)
    y = unsignedToSigned(i % 4096, 2048)
    i = int((i - y) / 4096)
    z = unsignedToSigned(i % 4096, 2048)
    return x,y,z

def limit(i, l, h):
    if(i > h):
        i = h
    if(i < l):
        i = l
    return i

def readU8(f):
    return ord(f.read(1))

def readU16(f):
    return ord(f.read(1))*256 + ord(f.read(1))

def readU32(f):
    return ord(f.read(1))*256*256*256 + ord(f.read(1))*256*256 + ord(f.read(1))*256 + ord(f.read(1))

def readS32(f):
    return unsignedToSigned(ord(f.read(1))*256*256*256 + ord(f.read(1))*256*256 + ord(f.read(1))*256 + ord(f.read(1)), 2**31)

usagetext = """minetestmapper.py [options]
  -i/--input <world_path>
  -o/--output <output_image.png>
  --bgcolor <color>
  --scalecolor <color>
  --playercolor <color>
  --origincolor <color>
  --drawscale
  --drawplayers
  --draworigin
  --drawunderground
Color format: '#000000'"""

def usage():
    print(usagetext)

try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:", ["help", "input=",
        "output=", "bgcolor=", "scalecolor=", "origincolor=",
        "playercolor=", "draworigin", "drawplayers", "drawscale",
        "drawunderground"])
except getopt.GetoptError as err:
    # print help information and exit:
    print(str(err))  # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

path = None
output = "map.png"
border = 0
scalecolor = "black"
bgcolor = ImageColor.getrgb("white")
origincolor = "red"
playercolor = "red"
drawscale = False
drawplayers = False
draworigin = False
drawunderground = False

sector_xmin = -1500 / 16
sector_xmax = 1500 / 16
sector_zmin = -1500 / 16
sector_zmax = 1500 / 16

for o, a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-i", "--input"):
        path = a
    elif o in ("-o", "--output"):
        output = a
    elif o == "--bgcolor":
        bgcolor = ImageColor.getrgb(a)
    elif o == "--scalecolor":
        scalecolor = ImageColor.getrgb(a)
    elif o == "--playercolor":
        playercolor = ImageColor.getrgb(a)
    elif o == "--origincolor":
        origincolor = ImageColor.getrgb(a)
    elif o == "--drawscale":
        drawscale = True
        border = 40
    elif o == "--drawplayers":
        drawplayers = True
    elif o == "--draworigin":
        draworigin = True
    elif o == "--drawunderground":
        drawunderground = True
    else:
        assert False, "unhandled option"

if path is None:
    print("Please select world path (eg. -i ../worlds/yourworld) (or use --help)")
    sys.exit(1)

if path[-1:] != "/" and path[-1:] != "\\":
    path = path + "/"

# Load color information for the blocks.
colors = {}
colors_str = {}
namenum=2
try:
    f = file("colors.txt")
except IOError:
    f = file(os.path.join(os.path.dirname(__file__), "colors.txt"))

for line in f:
    values = string.split(line)
    if len(values) < 4:
        continue
    identifier = values[0]
    is_hex = True
    for c in identifier:
        if c not in "0123456789abcdefABCDEF":
            is_hex = False
            break
    if is_hex:
        colors_str[int(values[0],16)] = namenum
        colors[namenum] = (
            int(values[1]),
            int(values[2]),
            int(values[3]))
    else:
        colors_str[values[0]] = namenum
        colors[namenum] = (
            int(values[1]),
            int(values[2]),
            int(values[3]))
    namenum+=1
f.close()

#print("colors: "+repr(colors))
#sys.exit(1)

xlist = []
zlist = []

# List all sectors to memory and calculate the width and heigth of the
# resulting picture.

conn = None
cur = None
if os.path.exists(path + "map.sqlite"):
    import sqlite3
    conn = sqlite3.connect(path + "map.sqlite")
    cur = conn.cursor()

    cur.execute("SELECT `pos` FROM `blocks`")
    while True:
        r = cur.fetchone()
        if not r:
            break

        x, y, z = getIntegerAsBlock(r[0])

        if x < sector_xmin or x > sector_xmax:
            continue
        if z < sector_zmin or z > sector_zmax:
            continue

        xlist.append(x)
        zlist.append(z)

if os.path.exists(path + "sectors2"):
    for filename in os.listdir(path + "sectors2"):
        for filename2 in os.listdir(path + "sectors2/" + filename):
            x = hex_to_int(filename)
            z = hex_to_int(filename2)
            if x < sector_xmin or x > sector_xmax:
                continue
            if z < sector_zmin or z > sector_zmax:
                continue
            xlist.append(x)
            zlist.append(z)

if os.path.exists(path + "sectors"):
    for filename in os.listdir(path + "sectors"):
        x = hex4_to_int(filename[:4])
        z = hex4_to_int(filename[-4:])
        if x < sector_xmin or x > sector_xmax:
            continue
        if z < sector_zmin or z > sector_zmax:
            continue
        xlist.append(x)
        zlist.append(z)

if len(xlist) == 0 or len(zlist) == 0:
    print("World does not exist.")
    sys.exit(1)

# Get rid of doubles
xlist, zlist = zip(*sorted(set(zip(xlist, zlist))))

minx = min(xlist)
minz = min(zlist)
maxx = max(xlist)
maxz = max(zlist)

w = (maxx - minx) * 16 + 16
h = (maxz - minz) * 16 + 16

print("Result image (w=" + str(w) + " h=" + str(h) + ") will be written to "
        + output)

im = Image.new("RGB", (w + border, h + border), bgcolor)
draw = ImageDraw.Draw(im)
#impix = im.load()

#stuff = {}
stuff_height = numpy.zeros([w,h],dtype = 'i2')
stuff_content = numpy.zeros([w,h],dtype='u2')
stuff_water = numpy.zeros([w,h],dtype = 'u2')
stuff_dnd = numpy.zeros([w,h],dtype=bool)
stuff = {'height':stuff_height,
    'content':stuff_content,
    'water':stuff_water,
    'dnd':stuff_dnd}

unknown_node_names = []
unknown_node_ids = []

starttime = time.time()

CONTENT_WATER = 2

def content_is_ignore(d):
    return d == 0
    #return d in [0, "ignore"]

def content_is_water(d):
    return (d == 2) | (d == 9)
    #return d in [2, 9]

def content_is_air(d):
    return (d == 126) | (d == 127) | (d == 254)
#    return d in [126, 127, 254, "air"]

def read_content(mapdata, version, datapos=None):
    if datapos==None:
        if version >= 24:
            mapdata = numpy.array(mapdata)
            x=numpy.arange(4096)
            return (mapdata[x*2] << 8) | (mapdata[x*2 + 1])

    if version >= 24:
        return (mapdata[datapos*2] << 8) | (mapdata[datapos*2 + 1])
    elif version >= 20:
        if mapdata[datapos] < 0x80:
            return mapdata[datapos]
        else:
            return (mapdata[datapos] << 4) | (mapdata[datapos + 0x2000] >> 4)
    elif 16 <= version < 20:
        return TRANSLATION_TABLE.get(mapdata[datapos], mapdata[datapos])
    else:
        raise Exception("Unsupported map format: " + str(version))


def read_mapdata(mapdata, version, plist, cdata, hdata, dnddata, day_night_differs, id_map, ignore, air):
    if(len(mapdata) < 4096):
        print("bad: " + xhex + "/" + zhex + "/" + yhex + " " + \
            str(len(mapdata)))
    else:
        chunkypos = ypos * 16
        mapdata = mapdata[:4096]
        mapdata = id_map[mapdata]
        if (mapdata==ignore).all():
            return plist
        mapdata = numpy.swapaxes(mapdata.reshape(16,16,16),1,0)
        mapdata = numpy.swapaxes(mapdata,1,2).reshape(16,256)
        y=15
        while len(plist)>0 and y>=0:
            content = mapdata[y][plist]
#            content = id_func(mapdata[y][plist])
#            watercontent = content_is_water(content)
#            wdata[plist] += watercontent
#            opaques = ~( (content_is_air(content) | content_is_ignore(content) | watercontent))
            opaques = ~( (content == ignore) | (content == air) )
            po = plist[opaques]
            cdata[po] = content[opaques]
            hdata[po] = chunkypos + y
            dnddata[po] = day_night_differs
            plist = plist[~opaques]
            y-=1
#        cdata[:] = id_map(cdata)
    return plist


read_map_time = 0

# Go through all sectors.
for n in range(len(xlist)):
    #if n > 500:
    #   break
    if n % 200 == 0:
        nowtime = time.time()
        dtime = nowtime - starttime
        try:
            n_per_second = 1.0 * n / dtime
        except ZeroDivisionError:
            n_per_second = 0
        if n_per_second != 0:
            seconds_per_n = 1.0 / n_per_second
            time_guess = seconds_per_n * len(xlist)
            remaining_s = time_guess - dtime
            remaining_minutes = int(remaining_s / 60)
            remaining_s -= remaining_minutes * 60
            print("Processing sector " + str(n) + " of " + str(len(xlist))
                    + " (" + str(round(100.0 * n / len(xlist), 1)) + "%)"
                    + " (ETA: " + str(remaining_minutes) + "m "
                    + str(int(remaining_s)) + "s)")

    xpos = xlist[n]
    zpos = zlist[n]

    xhex = int_to_hex3(xpos)
    zhex = int_to_hex3(zpos)
    xhex4 = int_to_hex4(xpos)
    zhex4 = int_to_hex4(zpos)

    sector1 = xhex4.lower() + zhex4.lower()
    sector2 = xhex.lower() + "/" + zhex.lower()

    ylist = []

    sectortype = ""

    if cur:
        psmin = getBlockAsInteger((xpos, -2048, zpos))
        psmax = getBlockAsInteger((xpos, 2047, zpos))
        cur.execute("SELECT `pos` FROM `blocks` WHERE `pos`>=? AND `pos`<=? AND (`pos` - ?) % 4096 = 0", (psmin, psmax, psmin))
        while True:
            r = cur.fetchone()
            if not r:
                break
            pos = getIntegerAsBlock(r[0])[1]
            ylist.append(pos)
            sectortype = "sqlite"
    try:
        for filename in os.listdir(path + "sectors/" + sector1):
            if(filename != "meta"):
                pos = int(filename, 16)
                if(pos > 32767):
                    pos -= 65536
                ylist.append(pos)
                sectortype = "old"
    except OSError:
        pass

    if sectortype == "":
        try:
            for filename in os.listdir(path + "sectors2/" + sector2):
                if(filename != "meta"):
                    pos = int(filename, 16)
                    if(pos > 32767):
                        pos -= 65536
                    ylist.append(pos)
                    sectortype = "new"
        except OSError:
            pass

    if sectortype == "":
        continue

    ylist.sort()

    # Make a list of pixels of the sector that are to be looked for.
    cdata = numpy.zeros(256,dtype='i4')
    hdata = numpy.zeros(256,dtype='i4')
    wdata = numpy.zeros(256,dtype='i4')
    dnddata = numpy.zeros(256,dtype=bool)
    plist = numpy.arange(256)

    # Go through the Y axis from top to bottom.
    for ypos in reversed(ylist):
        try:
            #print("("+str(xpos)+","+str(ypos)+","+str(zpos)+")")

            yhex = int_to_hex4(ypos)

            if sectortype == "sqlite":
                ps = getBlockAsInteger((xpos, ypos, zpos))
                cur.execute("SELECT `data` FROM `blocks` WHERE `pos`==? LIMIT 1", (ps,))
                r = cur.fetchone()
                if not r:
                    continue
                f = cStringIO.StringIO(r[0])
            else:
                if sectortype == "old":
                    filename = path + "sectors/" + sector1 + "/" + yhex.lower()
                else:
                    filename = path + "sectors2/" + sector2 + "/" + yhex.lower()
                f = file(filename, "rb")

            # Let's just memorize these even though it's not really necessary.
            version = readU8(f)
            flags = f.read(1)

            #print("version="+str(version))
            #print("flags="+str(version))

            # Check flags
            is_underground = ((ord(flags) & 1) != 0)
            day_night_differs = ((ord(flags) & 2) != 0)
            lighting_expired = ((ord(flags) & 4) != 0)
            generated = ((ord(flags) & 8) != 0)

            #print("is_underground="+str(is_underground))
            #print("day_night_differs="+str(day_night_differs))
            #print("lighting_expired="+str(lighting_expired))
            #print("generated="+str(generated))

            if version >= 22:
                content_width = readU8(f)
                params_width = readU8(f)

            # Node data
            dec_o = zlib.decompressobj()
            try:
                s = dec_o.decompress(f.read())
                mapdata = numpy.fromstring(s,">u2")
            except:
                mapdata = []
#                mapdata2 = numpy.array([])

            # Reuse the unused tail of the file
            f.close();
            f = cStringIO.StringIO(dec_o.unused_data)
            #print("unused data: "+repr(dec_o.unused_data))

            # zlib-compressed node metadata list
            dec_o = zlib.decompressobj()
            try:
                s=dec_o.decompress(f.read())
                metaliststr = numpy.fromstring(s,"u1")
                # And do nothing with it
            except:
                metaliststr = []

            # Reuse the unused tail of the file
            f.close();
            f = cStringIO.StringIO(dec_o.unused_data)
            #print("* dec_o.unused_data: "+repr(dec_o.unused_data))
            data_after_node_metadata = dec_o.unused_data

            if version <= 21:
                # mapblockobject_count
                readU16(f)

            if version == 23:
                readU8(f) # Unused node timer version (always 0)
            if version == 24:
                ver = readU8(f)
                if ver == 1:
                    num = readU16(f)
                    for i in range(0,num):
                        readU16(f)
                        readS32(f)
                        readS32(f)

            static_object_version = readU8(f)
            static_object_count = readU16(f)
            for i in range(0, static_object_count):
                # u8 type (object type-id)
                object_type = readU8(f)
                # s32 pos_x_nodes * 10000
                pos_x_nodes = readS32(f)/10000
                # s32 pos_y_nodes * 10000
                pos_y_nodes = readS32(f)/10000
                # s32 pos_z_nodes * 10000
                pos_z_nodes = readS32(f)/10000
                # u16 data_size
                data_size = readU16(f)
                # u8[data_size] data
                data = f.read(data_size)

            timestamp = readU32(f)
            #print("* timestamp="+str(timestamp))

            id_to_name = {}
            name_to_id = {}
            air = 1
            ignore = 0
            if version >= 22:
                name_id_mapping_version = readU8(f)
                num_name_id_mappings = readU16(f)
                #print("* num_name_id_mappings: "+str(num_name_id_mappings))
                for i in range(0, num_name_id_mappings):
                    node_id = readU16(f)
                    name_len = readU16(f)
                    name = f.read(name_len)
                    try:
                        id_to_name[node_id] = colors_str[name]
                    except:
                        ##TODO: Add to list of unknown colors
                        id_to_name[node_id] = 0
                    if name == 'air':
                        air = id_to_name[node_id]
                    if name == 'ignore':
                        ignore = id_to_name[node_id]
            if len(id_to_name)==0:
                id_map = [0,1]
            else:
                id_map = numpy.array([id_to_name[i] for i in sorted(id_to_name)])

            # Node timers
            if version >= 25:
                timer_size = readU8(f)
                num = readU16(f)
                for i in range(0,num):
                    readU16(f)
                    readS32(f)
                    readS32(f)
            plist = read_mapdata(mapdata, version, plist, cdata, hdata, dnddata, day_night_differs, id_map, ignore, air)
            # After finding all the pixels in the sector, we can move on to
            # the next sector without having to continue the Y axis.
            if len(plist) == 0 or y==ylist[0]:
                chunkxpos = (xpos-minx)*16
                chunkzpos = (zpos-minz)*16
                pos = (slice(chunkxpos,chunkxpos+16),slice(chunkzpos,chunkzpos+16))
                stuff['height'][pos] = hdata.reshape(16,16)
                stuff['content'][pos] = cdata.reshape(16,16)
                stuff['water'][pos] = wdata.reshape(16,16)
                stuff['dnd'][pos] = dnddata.reshape(16,16)
                break
        except Exception as e:
            print("Error at ("+str(xpos)+","+str(ypos)+","+str(zpos)+"): "+str(e))
            sys.stdout.write("Block data: ")
            for c in r[0]:
                sys.stdout.write("%2.2x "%ord(c))
            sys.stdout.write(os.linesep)
            sys.stdout.write("Data after node metadata: ")
            for c in data_after_node_metadata:
                sys.stdout.write("%2.2x "%ord(c))
            sys.stdout.write(os.linesep)
            traceback.print_exc()

print("Drawing image")
# Drawing the picture
starttime = time.time()
count_dnd=0
count_height=0
count_zero=0

c = stuff['content']
dnd = stuff['dnd']
hgh = stuff['height']
c0 = c[1:,:-1]
c1 = c[:-1,1:]
c2 = c[1:, 1:]
dnd0 = dnd[1:,:-1]
dnd1 = dnd[:-1,1:]
dnd2 = dnd[1:, 1:]
h0 = hgh[1:,:-1]
h1 = hgh[:-1,1:]
h2 = hgh[1:, 1:]
drop = (2*h0 - h1 - h2) * 12
drop = numpy.clip(drop,-255,32)

colors = numpy.array([bgcolor,bgcolor]+[colors[c] for c in sorted(colors)],dtype = 'i2')

pix = colors[stuff['content']]
pix[1:,:-1] += drop[:,:,numpy.newaxis]
pix = numpy.clip(pix,0,255)
pix = numpy.array(pix,dtype = 'u1')
impix = Image.fromarray(pix,'RGB')
impix = impix.transpose(Image.ROTATE_90)
im.paste(impix,(border,border))

if draworigin:
    draw.ellipse((minx * -16 - 5 + border, h - minz * -16 - 6 + border,
        minx * -16 + 5 + border, h - minz * -16 + 4 + border),
        outline=origincolor)

font = ImageFont.load_default()

if drawscale:
    draw.text((24, 0), "X", font=font, fill=scalecolor)
    draw.text((2, 24), "Z", font=font, fill=scalecolor)

    for n in range(int(minx / -4) * -4, maxx, 4):
        draw.text((minx * -16 + n * 16 + 2 + border, 0), str(n * 16),
            font=font, fill=scalecolor)
        draw.line((minx * -16 + n * 16 + border, 0,
            minx * -16 + n * 16 + border, border - 1), fill=scalecolor)

    for n in range(int(maxz / 4) * 4, minz, -4):
        draw.text((2, h - 1 - (n * 16 - minz * 16) + border), str(n * 16),
            font=font, fill=scalecolor)
        draw.line((0, h - 1 - (n * 16 - minz * 16) + border, border - 1,
            h - 1 - (n * 16 - minz * 16) + border), fill=scalecolor)

if drawplayers:
    try:
        for filename in os.listdir(path + "players"):
            f = file(path + "players/" + filename)
            lines = f.readlines()
            name = ""
            position = []
            for line in lines:
                p = string.split(line)
                if p[0] == "name":
                    name = p[2]
                    print(filename + ": name = " + name)
                if p[0] == "position":
                    position = string.split(p[2][1:-1], ",")
                    print(filename + ": position = " + p[2])
            if len(name) > 0 and len(position) == 3:
                x = (int(float(position[0]) / 10 - minx * 16))
                z = int(h - (float(position[2]) / 10 - minz * 16))
                draw.ellipse((x - 2 + border, z - 2 + border,
                    x + 2 + border, z + 2 + border), outline=playercolor)
                draw.text((x + 2 + border, z + 2 + border), name,
                    font=font, fill=playercolor)
            f.close()
    except OSError:
        pass

print("Saving")
im.save(output)

if unknown_node_names:
    sys.stdout.write("Unknown node names:")
    for name in unknown_node_names:
        sys.stdout.write(" "+name)
    sys.stdout.write(os.linesep)
if unknown_node_ids:
    sys.stdout.write("Unknown node ids:")
    for node_id in unknown_node_ids:
        sys.stdout.write(" "+str(hex(node_id)))
    sys.stdout.write(os.linesep)
