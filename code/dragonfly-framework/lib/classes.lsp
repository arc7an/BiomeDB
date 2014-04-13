;  Copyright (C) <2009> <Marc Hildmann, Greg Slepak>
;
;  This program is free software: you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation, either version 3 of the License, or
;  (at your option) any later version.
;
;  This program is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;  You should have received a copy of the GNU General Public License
;  along with this program.  If not, see <http://www.gnu.org/licenses/>.
;

(new Class 'Route)

(context Route)
(define (matches?) nil)
(define (run) nil)

(context 'Resource)
(define (Resource:Resource)
	(Dragonfly:die (context) " has no default action!")
)
(define (catch-all action id rformat) ; 'format' is a protected symbol
	(Dragonfly:die (context) ":" action " not defined!")
)

(context MAIN)
