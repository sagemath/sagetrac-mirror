libspatialindex
======================================

Description
-----------

The purpose of this library is to provide:
An extensible framework that will support robust spatial indexing methods.
Support for sophisticated spatial queries. Range, point location, nearest neighbor and k-nearest neighbor as well as parametric queries (defined by spatial constraints) should be easy to deploy and run.
Easy to use interfaces for inserting, deleting and updating information.
Wide variety of customization capabilities. Basic index and storage characteristics like the page size, node capacity, minimum fan-out, splitting algorithm, etc. should be easy to customize.
Index persistence. Internal memory and external memory structures should be supported. Clustered and non-clustered indices should be easy to be persisted.

Generic main memory and disk based storage managers.
R*-tree index (also supports linear and quadratic splitting).
MVR-tree index (a.k.a. PPR-tree).
TPR-tree index.
Advanced query capabilities, using Strategy and Visitor patterns.
Arbitrary shaped range queries, by defining generic geometry interfaces.
Large parameterization capabilities, including dimensionality, fill factor, node capacity, etc.
STR packing / bulk loading.

License
-------

libspatialindex changed from a LGPL to a MIT license as of the 1.8.0 release. For most situations, this should have no impact on the library’s use, but it should open it up for usage in situations that otherwise might have been problematic. Versions of libspatialindex prior to 1.8.0 were licensed LGPL 2.0, with the license description on this file. The codebase has been been updated, with licensing information replaced in headers and source files, to use the MIT license as of the 1.8.0+ release.
This change was made to support the inclusion of software depending on libspatialindex in static linking-only environments such as embedded systems and Apple’s iOS. libspatialindex versions prior to 1.8.0 will continue to live on as LGPL software, and developers can continue to contribute to them under terms of that license, but the main development effort, and ongoing maintenance, releases, and bug applications, will move forward using the new MIT license at http://github.com/libspatialindex/libspatialindex

Upstream Contact
----------------

https://libspatialindex.org/en/latest/