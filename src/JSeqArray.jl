# ===========================================================================
#
# JSeqArray.jl: Julia Interface to SeqArray Files
#
# Copyright (C) 2017    Xiuwen Zheng
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License Version 3 as
# published by the Free Software Foundation.
#
# JSeqArray is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with JSeqArray.
# If not, see <http://www.gnu.org/licenses/>.

module JSeqArray

using jugds

import Base: joinpath, show, print_with_color, println
import jugds: type_gdsfile

export	TypeSeqArray



####  Open and initialize the SeqArray binary library  ####

@static if is_apple()
	libfn = "libJSeqArray.dylib"
elseif is_windows()
	libfn = "libJSeqArray.dll"
elseif is_unix()
	libfn = "libJSeqArray.so"
else
	error("The platform is not supported.")
end

global libname = joinpath(Pkg.dir(), "JSeqArray", "deps", libfn)

if !isfile(libname)
	error("The SeqArray library cannot be found, please try Pkg.build(\"JSeqArray\").")
end

const LibSeqArray = libname

function __init__()
	ccall((:Init_GDS_Routines, LibSeqArray), Void, (Ptr{Void},), jugds.lib_c_api)
end



####  Type of GDS File and Node	 ####

type TypeSeqArray
	x :: type_gdsfile
end

end
