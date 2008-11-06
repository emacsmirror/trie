
;;; trie.el --- trie package


;; Copyright (C) 2004-2007 Toby Cubitt

;; Author: Toby Cubitt <toby-predictive@dr-qubit.org>
;; Version: 0.1
;; Keywords: trie, ternary search tree, completion
;; URL: http://www.dr-qubit.org/emacs.php


;; This file is NOT part of Emacs.
;;
;; This program is free software; you can redistribute it and/or
;; modify it under the terms of the GNU General Public License
;; as published by the Free Software Foundation; either version 2
;; of the License, or (at your option) any later version.
;;
;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program; if not, write to the Free Software
;; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
;; MA 02110-1301, USA.


;;; Commentary:
;;
;; Quick Overview
;; --------------
;; A trie is a data structure used to store keys that are ordered
;; sequences of elements (vectors, lists or strings in Elisp), in such a
;; way that both storage and retrieval are reasonably space- and
;; time-efficient. But, more importantly, searching for keys that match
;; various patterns can also be done efficiently. For example, returning
;; all strings with a given prefix, or searching for keys matching a
;; pattern containing wildcards, or searching for all keys within a given
;; Lewenstein distance of given string (though the latter two are not yet
;; implemented in this package - code contributions welcome!).
;;
;; You create a ternary search tree using `trie-create', create an
;; association using `trie-insert', retrieve an association using
;; `trie-lookup', find completions of a sequence using `trie-complete',
;; and map over a tree using `trie-map', `trie-mapc', `trie-mapcar', or
;; `trie-mapf'. Using `trie-stack', you can create an object that allows
;; the contents of the trie to be used like a stack; `trie-stack-pop'
;; pops elements off the stack one-by-one, whilst `trie-stack-push'
;; pushes things onto the stack.
;;
;; Note that there are two uses for a trie: as a lookup table, in which
;; case only the presence or absence of a key in the trie is significant,
;; or as an associative array, in which case each key carries some
;; associated data. Libraries for other data structure often only
;; implement lookup tables, leaving it up to you to implement an
;; associative array on top of this (by storing key+data pairs in the
;; data structure's keys, then defining a comparison function that only
;; compares the key part). However, for a trie, this would be slightly
;; less space-efficient than it needs to be, so this package does the
;; opposite: it implements associative arrays, and leaves it up to you to
;; use them as lookup tables if you so desire (with no loss of
;; space-efficiency).
;;
;;
;; Different Types of Trie
;; -----------------------
;; There are numerous ways to implement trie data structures internally,
;; each with its own trade-offs. By viewing a trie as a tree whose nodes
;; are themselves lookup tables for key elements, this package is able to
;; support all types of trie, providing there exists (or you write!) an
;; Elisp implementation of the corresponding type of lookup table. The
;; best implementation will depend on what trade-offs are appropriate for
;; your particular application. The following gives an overview of the
;; advantages and disadvantages of various types of trie. (Not all of the
;; underlying lookup tables have been implemented in Elisp yet, so using
;; some of them would require writing the missing Elisp package!)
;;
;; One of the most effective all-round implementations of a trie is a
;; ternary search tree, which can be viewed as a tree of binary trees. If
;; basic binary search trees are used for the nodes of the trie, we get a
;; basic ternary search tree. If self-balancing binary trees are used
;; (e.g. AVL or red-black trees), we get a self-balancing ternary search
;; tree. If splay trees are used, we get yet another self-organising
;; variant of a ternary search tree. All ternary search trees have, in
;; common, good space-efficiency. The time-efficiencies for the various
;; trie operations are also good, assuming the underlying binary trees
;; are balanced. Under that assumption, all variants of ternary search
;; trees described below have the same asymptotic time-complexity for all
;; trie operations.
;;
;; Self-balancing trees ensure the underlying binary trees are always
;; close to perfectly balanced, with the usual trade-offs between the
;; different the types of self-balancing binary tree: AVL trees are
;; slightly more efficient for lookup operations than red-black trees,
;; but are slightly less efficienct for insertion operations, and less
;; efficient for deletion operations. Splay trees give good average-case
;; complexity and are simpler to implement than AVL or red-black trees
;; (which can mean they're faster in practice!), at the expense of poor
;; worst-case complexity.
;;
;; If your tries are going to be static (i.e. created once and rarely
;; modified), then using perfectly balanced binary search trees might be
;; more appropriate. Perfectly balancing the binary trees is very
;; inefficient, but it only has to be when the trie is first created or
;; modified. Lookup operations will then be as efficient as possible for
;; ternary search trees, and the implementation will be much simpler (so
;; probably faster) than a self-balancing tree, without the space and
;; time overhead required to keep track of rebalancing.
;;
;; On the other hand, adding data to a binary search tree in a random
;; order usually results in a reasonably balanced tree. If this is the
;; likely scenario, using a basic binary tree without bothering to
;; balance it at all might be quite efficient, and, being even simpler to
;; implement, could be faster overall.
;;
;; A digital trie is a different implementation of a trie, which can be
;; viewed as a tree of arrays, and has different space- and
;; time-complexities than a ternary search tree. Roughly speaking, a
;; digital trie has worse space-complexity, but better
;; time-complexity. Using hash tables instead of arrays for the nodes
;; gives something similar to a digital trie, potentially with better
;; space-complexity and the same amortised time-complexity, but at the
;; expense of occasional significant inefficiency when inserting and
;; deleting (whenever the hash table has to be resized). Indeed, an array
;; can be viewed as a perfect hash table, but as such it requires the
;; number of possible values to be known in advance.
;;
;; Finally, if you really need optimal efficiency from your trie, you
;; could even write a custom type of underlying lookup table, optimised
;; for your specific needs.
;;
;;
;; This package uses the AVL tree package avl-tree.el and the heap
;; package heap.el.


;;; Change Log:
;;
;; Version 0.1
;; * Initial release (complete rewrite from scratch of tstree.el!)
;; * Ternary search trees are now implemented as a tree of avl trees, which
;;   has numerous advantages: self-balancing trees guarantee O(log n)
;;   complexity regardless of how the tree is built; deletion is now done
;;   properly.
;; * unlike tstree.el, trie.el is general enough to implement all sorts
;;   of tries, not just ternary search trees (though these remain the
;;   default).
;; * Up to "tstree"->"trie" renaming, many functions are drop-in
;;   replacements for tstree.el functions. However, insertion and rank
;;   functions are no longer stored in the data structure, so
;;   corresponidng arguments are no longer optional. A single
;;   `trie-complete' function now deals with sorting completions in both
;;   lexical or arbitrary order, the ranking function being passed as an
;;   optional argument in the latter case. And functions can no longer
;;   operate over multiple data structures at once; i.e. they no longer
;;   accept lists of trees as arguments. (These features belong in higher
;;   level libraries, and the efficiency loss is negligible.)



;;; Code:

(eval-when-compile (require 'cl))
(require 'avl-tree)
(require 'heap)



;;; ================================================================
;;;                Replacements for CL functions

;; copied from cl-extra.el
(defun trie--subseq (seq start &optional end)
  "Return the subsequence of SEQ from START to END.
If END is omitted, it defaults to the length of the sequence.
If START or END is negative, it counts from the end."
  (if (stringp seq) (substring seq start end)
    (let (len)
      (and end (< end 0) (setq end (+ end (setq len (length seq)))))
      (when (< start 0)
	(setq start (+ start (or len (setq len (length seq))))))
      (cond ((listp seq)
	     (if (> start 0) (setq seq (nthcdr start seq)))
	     (if end
		 (let ((res nil))
		   (while (>= (setq end (1- end)) start)
		     (push (pop seq) res))
		   (nreverse res))
	       (copy-sequence seq)))
	    (t
	     (or end (setq end (or len (length seq))))
	     (let ((res (make-vector (max (- end start) 0) nil))
		   (i 0))
	       (while (< start end)
		 (aset res i (aref seq start))
		 (setq i (1+ i) start (1+ start)))
	       res))))))


(defsubst trie--seq-append (seq el)
  "Append EL to the end of sequence SEQ."
  (cond
   ((stringp seq) (concat seq (string el)))
   ((vectorp seq) (vconcat seq (vector el)))
   ((listp seq)	  (append seq (list el)))))



;;; ================================================================
;;;     Internal functions only for use within the trie package


;;; ----------------------------------------------------------------
;;;           Functions and macros for handling a trie.

;; symbol used to denote a trie leaf node
(defconst trie--terminator '--trie--terminator)

(defstruct
  (trie-
   :named
   (:constructor nil)
   (:constructor trie--create
		 (comparison-function &optional (type 'avl)
		  &aux
		  (createfun
		   (cond
		    ((eq type 'avl) 'avl-tree-create)
		    (t (error "trie--create: unknown trie TYPE, %s" type))))
		  (insertfun
		   (cond
		    ((eq type 'avl) 'avl-tree-enter)
		    (t (error "trie--create: unknown trie TYPE, %s" type))))
		  (deletefun
		   (cond
		    ((eq type 'avl) 'avl-tree-delete)
		    (t (error "trie--create: unknown trie TYPE, %s" type))))
		  (lookupfun
		   (cond
		    ((eq type 'avl) 'avl-tree-member)
		    (t (error "trie--create: unknown trie TYPE, %s" type))))
		  (mapfun
		   (cond
		    ((eq type 'avl) 'avl-tree-mapc)
		    (t (error "trie--create: unknown trie TYPE, %s" type))))
		  (emptyfun
		   (cond
		    ((eq type 'avl) 'avl-tree-empty)
		    (t (error "trie--create: unknown trie TYPE, %s" type))))
		  (stack-createfun
		   (cond
		    ((eq type 'avl) 'avl-tree-stack)
		    (t (error "trie--create: unknown trie TYPE, %s" type))))
		  (stack-popfun
		   (cond
		    ((eq type 'avl) 'avl-tree-stack-pop)
		    (t (error "trie--create: unknown trie TYPE, %s" type))))
		  (stack-emptyfun
		   (cond
		    ((eq type 'avl) 'avl-tree-stack-empty-p)
		    (t (error "trie--create: unknown trie TYPE, %s" type))))
		  (transform-for-print
		   (cond
		    ((eq type 'avl) 'trie--avl-transform-for-print)
		    (t (error "trie--create: unknown trie TYPE, %s" type))))
		  (transform-from-read
		   (cond
		    ((eq type 'avl) 'trie--avl-transform-from-read)
		    (t (error "trie--create: unknown trie TYPE, %s" type))))
		  (cmpfun (trie--wrap-cmpfun comparison-function))
		  (root (trie--node-create-root createfun cmpfun))
		  ))
   (:constructor trie--create-custom
		 (comparison-function
		  &key
		  (createfun 'avl-tree-create-bare)
		  (insertfun 'avl-tree-enter)
		  (deletefun 'avl-tree-delete)
		  (lookupfun 'avl-tree-member)
		  (mapfun 'avl-tree-mapc)
		  (emptyfun 'avl-tree-empty)
		  (stack-createfun 'avl-tree-stack)
		  (stack-popfun 'avl-tree-stack-pop)
		  (stack-emptyfun 'avl-tree-stack-empty-p)
		  (transform-for-print 'trie--avl-transform-for-print)
		  (transform-from-read 'trie--avl-transform-from-read)
		  &aux
		  (cmpfun (trie--wrap-cmpfun comparison-function))
		  (root (trie--node-create-root createfun cmpfun))
		  ))
   (:copier nil))
  root comparison-function cmpfun
  createfun insertfun deletefun lookupfun mapfun emptyfun
  stack-createfun stack-popfun stack-emptyfun
  transform-for-print transform-from-read print-form)


(defun trie--wrap-cmpfun (cmpfun)
  ;; wrap CMPFUN for use in a subtree
  `(lambda (a b)
     (setq a (trie--node-split a)
	   b (trie--node-split b))
     (cond ((eq a trie--terminator)
	    (if (eq b trie--terminator) nil t))
	   ((eq b trie--terminator) nil)
	   (t (,cmpfun a b)))))


;;; ----------------------------------------------------------------
;;;          Functions and macros for handling a trie node.

(defstruct
  (trie--node
   (:type vector)
   (:constructor nil)
   (:constructor trie--node-create
		 (split trie
		  &aux (subtree (funcall (trie--createfun trie)
					 (trie--cmpfun trie)))))
   (:constructor trie--node-create-data
		 (data &aux (split trie--terminator) (subtree data)))
   (:constructor trie--node-create-dummy
		 (split &aux (subtree nil)))
   (:constructor trie--node-create-root
		 (createfun cmpfun
		  &aux
		  (split nil)
		  (subtree (let ((trie-depth 0))
			     (funcall createfun cmpfun)))))
   (:copier nil))
   split subtree)

;; data is stored in the subtree cell of a terminal node
(defalias 'trie--node-data 'trie--node-subtree)

(defsetf trie--node-data (node) `(setf (trie--node-subtree ,node)))

(defmacro trie--node-data-p (node)
  ;; Return t if NODE is a data node, nil otherwise.
  `(eq (trie--node-split ,node) trie--terminator))

(defmacro trie--node-p (node)
  ;; Return t if NODE is a trie--node, nil otherwise.
  ;; Have to define this ourselves, because we created a defstruct
  ;; without any identifying tags (i.e. (:type vector)) for efficiency.
  `(and (vectorp ,node)
	(= (length ,node) 2)
	(or (trie--node-data-p ,node)
	    (trie--p (trie--node-subtree ,node)))))


(defun trie--node-find (trie sequence)
  ;; Returns the node corresponding to SEQUENCE, or nil if none found.
  (let ((node (trie--root trie))
	(len (length sequence))
	(i -1))
    ;; descend trie until we find SEQUENCE or run out of trie
    (while (and node (< (incf i) len))
      (setq node
	    (funcall (trie--lookupfun trie)
		     (trie--node-subtree node)
		     (trie--node-create-dummy (elt sequence i))
		     nil)))
    node))


(defmacro trie--find-data-node (node trie)
  ;; Return data node from NODE's subtree, or nil if NODE has no data node in
  ;; its subtree.
  `(funcall (trie--lookupfun ,trie)
	    (trie--node-subtree ,node)
	    (trie--node-create-dummy trie--terminator)
	    nil))


(defmacro trie--find-data (node trie)
  ;; Return data associated with sequence corresponding to NODE, or nil if
  ;; sequence has no associated data.
  `(let ((node (trie--find-data-node ,node ,trie)))
     (when node (trie--node-data node))))



;;; ----------------------------------------------------------------
;;;          Functions and macros for handling trie-stacks

(defstruct (trie--stack
	    (:constructor nil)
	    (:constructor
	     trie--stack-create
	     (trie
	      &optional
	      (type 'vector)
	      reverse
	      &aux
	      (stack-createfun (trie--stack-createfun trie))
	      (stack-popfun (trie--stack-popfun trie))
	      (stack-emptyfun (trie--stack-emptyfun trie))
	      (store
	       (if (trie-empty trie)
		   nil
		 (list (cons
			(cond ((eq type 'list) ())
			      ((eq type 'string) "")
			      (t []))
			(funcall stack-createfun
				 (trie--node-subtree (trie--root trie))
				 reverse)))))
	      (pushed '())
	      ))
	    (:constructor
	     trie--completion-stack-create
	     (trie prefix
	      &optional
	      reverse
	      &aux
	      (stack-createfun (trie--stack-createfun trie))
	      (stack-popfun (trie--stack-popfun trie))
	      (stack-emptyfun (trie--stack-emptyfun trie))
	      (store (trie--completion-stack-construct-store
		      trie prefix reverse))
	      (pushed '())
	      ))
	    (:copier nil))
  reverse predicatefun stack-createfun stack-popfun stack-emptyfun
  store pushed)


(defun trie--completion-stack-construct-store (trie prefix reverse)
  ;; Construct store for completion stack based on TRIE.
  (let (accumulate node)
    (if (or (atom prefix)
	    (and (listp prefix)
		 (not (sequencep (car prefix)))))
	(setq prefix (list prefix))
      (setq prefix
	    (sort prefix
		  (trie-construct-sortfun
		   (trie--comparison-function trie)
		   (not reverse)))))
    (dolist (pfx prefix)
      (when (setq node (trie--node-find trie pfx))
	(push (cons pfx (funcall (trie--stack-createfun trie)
				 (trie--node-subtree node)
				 reverse))
	      accumulate)))
    accumulate))


(defun trie--stack-repopulate (stack)
  ;; Recursively push children of the node at the head of STACK onto the front
  ;; of STACK, until a data node is reached.

  ;; nothing to do if stack is empty
  (unless (trie-stack-empty-p stack)
    (let ((node (funcall (trie--stack-stack-popfun stack)
			 (cdar (trie--stack-store stack))))
	  (seq (caar (trie--stack-store stack))))
      (when (funcall (trie--stack-stack-emptyfun stack)
		     (cdar (trie--stack-store stack)))
	;; effectively (pop (trie--stack-store stack)) w/o compilter warnings
	(setf (trie--stack-store stack) (cdr (trie--stack-store stack))))

      (while (not (trie--node-data-p node))
	(push
	 (cons (trie--seq-append seq (trie--node-split node))
	       (funcall (trie--stack-stack-createfun stack)
			(trie--node-subtree node)))
	 (trie--stack-store stack))
	(setq node (funcall (trie--stack-stack-popfun stack)
			    (cdar (trie--stack-store stack)))
	      seq (caar (trie--stack-store stack)))
	(when (funcall (trie--stack-stack-emptyfun stack)
		       (cdar (trie--stack-store stack)))
	  ;; effectively (pop (trie--stack-store stack)) w/o compiler warnings
	  (setf (trie--stack-store stack) (cdr (trie--stack-store stack)))))

      (push (cons seq (trie--node-data node)) (trie--stack-store stack)))))



;;; ----------------------------------------------------------------
;;;              print/read transformation functions

(defun trie-transform-for-print (trie)
  "Transform TRIE to print form."
  (when (trie--transform-for-print trie)
    (if (trie--print-form trie)
	(error "Trie has already been transformed to print-form")
      (setf (trie--print-form trie) t)
      (funcall (trie--transform-for-print trie) trie))))


(defun trie-transform-from-read (trie)
  "Transform TRIE from print form."
  (when (trie--transform-from-read trie)
    (if (not (trie--print-form trie))
	(error "Trie is not in print-form")
      (setf (trie--print-form trie) nil)
      (funcall (trie--transform-from-read trie) trie))))


(defun trie--avl-transform-for-print (trie)
  ;; transform avl-tree based TRIE to print form.
  (trie-mapc-internal
   (lambda (avl seq) (setf (avl-tree--cmpfun avl) nil))
   trie))


(defun trie--avl-transform-from-read (trie)
  ;; transform avl-tree based TRIE from print form."
  (let ((--trie-avl-transform--cmpfun (trie--cmpfun trie)))
    (trie-mapc-internal
     (lambda (avl seq)
       (setf (avl-tree--cmpfun avl) --trie-avl-transform--cmpfun))
     trie)))



;;; ----------------------------------------------------------------
;;;               Miscelaneous internal macros

(defun trie--mapc (--trie--mapc--function --trie--mapc--mapfun
		   --trie--mapc--root --trie--mapc--seq
		   &optional --trie--mapc--reverse)
  ;; Apply TRIE--MAPC--FUNCTION to all elements in a trie beneath
  ;; TRIE--MAPC--ROOT, which should correspond to the sequence
  ;; TRIE--MAPC--SEQ. TRIE--MAPC--FUNCTION is passed two arguments: the trie
  ;; node itself and the sequence it corresponds to. It is applied in
  ;; ascending order, or descending order if TRIE--MAPC--REVERSE is non-nil.

  ;; The absurd argument names are to lessen the likelihood of dynamical
  ;; scoping bugs caused by a supplied function binding a variable with the
  ;; same name as one of the arguments.
  (funcall
   --trie--mapc--mapfun
   (lambda (--trie--mapc--node)
     ;; data node: apply function
     (if (trie--node-data-p --trie--mapc--node)
	 (funcall --trie--mapc--function --trie--mapc--node --trie--mapc--seq)
       ;; internal node: append split value to seq and keep descending
       (trie--mapc --trie--mapc--function --trie--mapc--mapfun
		   --trie--mapc--node
		   (trie--seq-append (copy-sequence --trie--mapc--seq)
				     (trie--node-split --trie--mapc--node))
		   --trie--mapc--reverse)))
   ;; --TRIE--MAPC--MAPFUN target
   (trie--node-subtree --trie--mapc--root)
   --trie--mapc--reverse))


(defun trie-mapc-internal (function trie &optional type)
  "Apply FUNCTION to all internal associative arrays within TRIE.
FUNCTION is passed two arguments: an associative array, and the
sequence it corresponds to.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence passed to function. Defaults to
vector."
  (trie--mapc-internal function (trie--mapfun trie) (trie--root trie)
		       (cond ((eq type 'string) "")
			     ((eq type 'lisp) ())
			     (t []))))


(defun trie--mapc-internal (--trie--mapc-internal--function
			     --trie--mapc-internal--mapfun
			     --trie--mapc-internal--root
			     --trie--mapc-internal--seq)
  (funcall
   --trie--mapc-internal--mapfun
   (lambda (--trie--mapc-internal--node)
     ;; data node
     (unless (trie--node-data-p --trie--mapc-internal--node)
       (funcall --trie--mapc-internal--function
		(trie--node-subtree --trie--mapc-internal--node)
		--trie--mapc-internal--seq)
       (trie--mapc-internal
	--trie--mapc-internal--function
	--trie--mapc-internal--mapfun
	--trie--mapc-internal--node
	(trie--seq-append (copy-sequence --trie--mapc-internal--seq)
			  (trie--node-split --trie--mapc-internal--node)))))
   (trie--node-subtree --trie--mapc-internal--root)))


(defmacro trie--complete-construct-accumulator (maxnum filter)
  ;; Does what it says on the tin! | sed -e 's/on/in/' -e 's/tin/macro name/'
  `(cond
    ((and ,filter ,maxnum)
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (when (funcall ,filter seq data)
	   (aset trie--complete-accumulate 0
		 (cons (cons seq data)
		       (aref trie--complete-accumulate 0)))
	   (and (>= (length (aref trie--complete-accumulate 0)) ,maxnum)
		(throw 'trie-complete--done nil))))))
    ((and (not ,filter) ,maxnum)
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (aset trie--complete-accumulate 0
	       (cons (cons seq data)
		     (aref trie--complete-accumulate 0)))
	 (and (>= (length (aref trie--complete-accumulate 0)) ,maxnum)
	      (throw 'trie-complete--done nil)))))
    ((and ,filter (not ,maxnum))
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (when (funcall ,filter seq data)
	   (aset trie--complete-accumulate 0
		 (cons (cons seq data)
		       (aref trie--complete-accumulate 0)))))))
    ((and (not ,filter) (not ,maxnum))
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (aset trie--complete-accumulate 0
	       (cons (cons seq data)
		     (aref trie--complete-accumulate 0))))))))


(defmacro trie--complete-construct-ranked-accumulator (maxnum filter)
  ;; Does what it says on the tin! | sed -e 's/on/in/' -e 's/tin/macro name/'
  `(cond
    ((and ,filter ,maxnum)
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (when (funcall ,filter seq data)
	   (heap-add trie--complete-accumulate (cons seq data))
	   (and (> (heap-size trie--complete-accumulate) ,maxnum)
		(heap-delete-root trie--complete-accumulate))))))
    ((and ,filter (not ,maxnum))
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (when (funcall ,filter seq data)
	   (heap-add trie--complete-accumulate (cons seq data))))))
    ((and (not ,filter) ,maxnum)
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (heap-add trie--complete-accumulate (cons seq data))
	 (and (> (heap-size trie--complete-accumulate) ,maxnum)
	      (heap-delete-root trie--complete-accumulate)))))
    ((and (not ,filter) (not ,maxnum))
     (lambda (node seq)
       (let ((data (trie--node-data node)))
	 (heap-add trie--complete-accumulate (cons seq data)))))))




;;; ================================================================
;;;        The public functions which operate on tries.

(defalias 'trie-create 'trie--create
  "Return a new trie that uses comparison function COMPARISON-FUNCTION.

A trie stores sequences (strings, vectors or lists) along with
associated data. COMPARISON-FUNCTEION should accept two
arguments, each being an element of such a sequence, and return t
if the first is strictly smaller than the second.

The optional argument TYPE specifies the type of trie to
create. However, the only one that is currently implemented is
the default, so this argument is useless. (See also
`trie-create-custom'.)")



(defalias 'trie-create-custom 'trie--create-custom
  "Return a new trie that uses comparison function COMPARISON-FUNCTION.

A trie stores sequences (strings, vectors or lists) along with
associated data. COMPARISON-FUNCTION should accept two arguments,
each being an element of such a sequence, and return t if the
first is strictly smaller than the second.

The remaining arguments: CREATEFUN, INSERTFUN, DELETEFUN,
LOOKUPFUN, MAPFUN, EMPTYFUN, STACK-CREATEFUN, STACK-POPFUN and
STACK-EMPTYFUN, determine the type of trie that is created.

CREATEFUN is called as follows:

  (CREATEFUN COMPARISON-FUNCTION)

and should return a data structure (\"ARRAY\") that can be used
as an associative array, where two elements A and B are equal if
the following is non-nil:

  (and (COMPARISON-FUNCTION b a)
       (COMPARISON-FUNCTION b a))

When CREATEFUN is called, the depth in the trie at which the
associative array is being created can be accessed via the
variable `trie-depth'. This can be used, for example, to create
hybrid trie structures in which different types of associative
array are used at different depths in the trie. (Note that, in
this case, all the other functions described below must be able
to correctly handle *any* of these types of associate array.)

INSERTFUN, DELETEFUN, LOOKUPFUN, MAPFUN and EMPTYFUN should
insert, delete, lookup, map over, and check-if-there-exist-any
elements in the associative array. They are called as follows:

  (INSERTFUN array element &optional updatefun)
  (DELETEFUN array element)
  (LOOKUPFUN array element &optional nilflag)
  (MAPFUN function array &optional reverse)
  (EMPTYFUN array)

INSERTFUN should return the new element, which will be ELEMENT
itself unless UPDATEFUN is specified. In the latter case, it
should pass two arguments to UPDATEFUN, ELEMENT and the matching
element in the associate array, and replace that element with the
return value.

LOOKUPFUN should return the element from the associative array
that is equal to ELEMENT, or NILFLAG if no match exists.

MAPFUN should map FUNCTION over all elements in the order defined by
COMPARISON-FUNCTION, or in reverse order if REVERSE is non-nil.


STACK-CREATEFUN, STACK-POPFUN and STACK-EMPTYFUN should allow the
associative array to be used as a stack. STACK-CREATEFUN is
called as follows:

  (STACK-CREATEFUN array)

and should return a data structure (\"STACK\") that behaves like
a sorted stack of all elements in the associative array. I.e.
successive calls to

  (STACK-POPFUN stack)

should return elements from the associative array in the order
defined by COMPARISON-FUNCTION.

  (STACK-EMPTYFUN stack)

should return non-nil if the stack is empty, nil otherwise.

Note: to avoid nasty dynamic scoping bugs, the supplied functions
must *not* bind any variables with names commencing \"--\".")



(defalias 'trie-comparison-function 'trie--comparison-function
  "Return the comparison function for TRIE.")


(defalias 'trie-p 'trie--p
  "Return t if argument is a trie, nil otherwise.")


(defun trie-empty (trie)
  "Return t if the TRIE is empty, nil otherwise."
  (if (trie--print-form trie)
      (error "Attempted to operate on trie that is in print-form")
    (funcall (trie--emptyfun trie)
	     (trie--node-subtree (trie--root trie)))))


(defun trie-construct-sortfun (cmpfun &optional reverse)
  "Construct function to compare key sequences, based on a CMPFUN
that compares individual elements of the sequence. Order is
reversed if REVERSE is non-nil."
  (if reverse
      `(lambda (a b)
	 (let (cmp)
	   (catch 'compared
	     (dotimes (i (min (length a) (length b)))
	       (cond ((,cmpfun (elt b i) (elt a i)) (throw 'compared t))
		     ((,cmpfun (elt a i) (elt b i)) (throw 'compared nil))))
	     (< (length a) (length b)))))
    `(lambda (a b)
       (let (cmp)
	 (catch 'compared
	   (dotimes (i (min (length a) (length b)))
	     (cond ((,cmpfun (elt a i) (elt b i)) (throw 'compared t))
		   ((,cmpfun (elt b i) (elt a i)) (throw 'compared nil))))
	   (< (length a) (length b)))))))



;;; ----------------------------------------------------------------
;;;                      Mapping over tries

(defun trie-map (function trie &optional type reverse)
  "Modify all elements in TRIE by applying FUNCTION to them.

FUNCTION should take two arguments: a sequence stored in the trie
and its associated data. Its return value replaces the existing
data.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence passed to FUNCTION. Defaults to
vector.

FUNCTION is applied in ascending order, or descending order if
REVERSE is non-nil.

Note: to avoid nasty dynamic scoping bugs, FUNCTION must *not*
bind any variables with names commencing \"--\"."
  (if (trie--print-form trie)
      (error "Attempted to operate on trie that is in print-form")
    (let ((--trie-map--function function)) ; try to avoid dynamic scoping bugs
      (trie--mapc
       (lambda (node seq)
	 (setf (trie--node-data node)
	       (funcall --trie-map--function seq (trie--node-data node))))
       (trie--mapfun trie)
       (trie--root trie)
       (cond ((eq type 'string) "") ((eq type 'lisp) ()) (t []))
       reverse))))


(defun trie-mapc (function trie &optional type reverse)
  "Apply FUNCTION to all elements in TRIE for side effect only.

FUNCTION should take two arguments: a sequence stored in the trie
and its associated data.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence passed to FUNCTION. Defaults to
vector.

FUNCTION is applied in ascending order, or descending order if
REVERSE is non-nil.

Note: to avoid nasty dynamic scoping bugs, FUNCTION must *not*
bind any variables with names commencing \"--\"."
  (if (trie--print-form trie)
      (error "Attempted to operate on trie that is in print-form")
    (let ((--trie-mapc--function function)) ; try to avoid dynamic scoping bugs
      (trie--mapc
       (lambda (node seq)
	 (funcall --trie-mapc--function seq (trie--node-data node)))
       (trie--mapfun trie)
       (trie--root trie)
       (cond ((eq type 'string) "") ((eq type 'lisp) ()) (t []))
       reverse))))


(defun trie-mapf (function combinator trie &optional type reverse)
  "Apply FUNCTION to all elements in TRIE, and combine the results
using COMBINATOR.

FUNCTION should take two arguments: a sequence stored in the
trie, and its associated data.

Optional argument TYPE (one of the symbols vector, lisp or
string; defaults to vector) sets the type of sequence passed to
FUNCTION. If TYPE is 'string, it must be possible to apply the
function `string' to the individual elements of key sequences
stored in TRIE.

The FUNCTION is applied and the results combined in ascending
order, or descending order if REVERSE is non-nil.

Note: to avoid nasty dynamic scoping bugs, FUNCTION and
COMBINATOR must *not* bind any variables with names
commencing \"--\"."
  (if (trie--print-form trie)
      (error "Attempted to operate on trie that is in print-form")
    (let ((--trie-mapf--function function) ; try to avoid dynamic scoping bugs
	  --trie-mapf--accumulate)
      (trie--mapc
       (lambda (node seq)
	 (setq --trie-mapf--accumulate
	       (funcall combinator
			(funcall --trie-mapf--function
				 seq (trie--node-data node))
			--trie-mapf--accumulate)))
       (trie--mapfun trie)
       (trie--root trie)
       (cond ((eq type 'string) "") ((eq type 'lisp) ()) (t []))
       reverse)
      --trie-mapf--accumulate)))


(defun trie-mapcar (function trie &optional type reverse)
  "Apply FUNCTION to all elements in TRIE, and make a list of the results.

FUNCTION should take two arguments: a sequence stored in the trie
and its associated data.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence passed to FUNCTION. Defaults to
vector.

The FUNCTION is applied and the list constructed in ascending
order, or descending order if REVERSE is non-nil.

Note that if you don't care about the order in which FUNCTION is
applied, just that the resulting list is in the correct order,
then

  (trie-mapf function 'cons trie type (not reverse))

is more efficient.

Note: to avoid nasty dynamic scoping bugs, FUNCTION must *not*
bind any variables with names commencing \"--\"."
  (if (trie--print-form trie)
      (error "Attempted to operate on trie that is in print-form")
    (nreverse (trie-mapf function 'cons trie type reverse))))



;;; ----------------------------------------------------------------
;;;                    Using tries as stacks

(defun trie-stack (trie &optional type reverse)
  "Return an object that allows TRIE to be accessed as if it were a stack.

The stack is sorted in \"lexical\" order, i.e. the order defined
by the trie's comparison function, or in reverse order if REVERSE
is non-nil. Calling `trie-stack-pop' pops the top element (a key
and its associated data) from the stack.

Optional argument TYPE (one of the symbols vector, lisp or
string) sets the type of sequence used for the keys.

Note that any modification to TRIE *immediately* invalidates all
trie-stacks created before the modification (in particular,
calling `trie-stack-pop' will give unpredictable results).

Operations on trie-stacks are significantly more efficient than
constructing a real stack from the trie and using standard stack
functions. As such, they can be useful in implementing efficient
algorithms on tries. However, in cases where mapping functions
`trie-mapc', `trie-mapcar' or `trie-mapf' would be sufficient, it
is better to use one of those instead."
  (if (trie--print-form trie)
      (error "Attempted to operate on trie that is in print-form")
    (let ((stack (trie--stack-create trie type reverse)))
      (trie--stack-repopulate stack)
      stack)))


(defun trie-complete-stack (trie prefix &optional reverse)
  "Return an object that allows completions of PREFIX to be accessed
as if they were a stack.

The stack is sorted in \"lexical\" order, i.e. the order defined
by TRIE's comparison function, or in reverse order if REVERSE is
non-nil. Calling `trie-stack-pop' pops the top element (a key and
its associated data) from the stack.

PREFIX must be a sequence (vector, list or string) that forms the
initial part of a TRIE key. (If PREFIX is a string, it must be
possible to apply `string' to individual elements of TRIE keys.)
The completions returned in the alist will be sequences of the
same type as KEY. If PREFIX is a list of sequences, completions
of all sequences in the list are included in the stack. All
sequences in the list must be of the same type.

Note that any modification to TRIE *immediately* invalidates all
trie-stacks created before the modification (in particular,
calling `trie-stack-pop' will give unpredictable results).

Operations on trie-stacks are significantly more efficient than
constructing a real stack from completions of PREFIX in TRIE and
using standard stack functions. As such, they can be useful in
implementing efficient algorithms on tries. However, in cases
where `trie-complete' or `trie-complete-ordered' is sufficient,
it is better to use one of those instead."
  (if (trie--print-form trie)
      (error "Attempted to operate on trie that is in print-form")
    (let ((stack (trie--completion-stack-create trie prefix reverse)))
      (trie--stack-repopulate stack)
      stack)))


(defun trie-stack-pop (trie-stack)
  "Pop the first element from TRIE-STACK.
Returns nil if the stack is empty."
  ;; if elements have been pushed onto the stack, pop those first
  (if (trie--stack-pushed trie-stack)
      (pop (trie--stack-pushed trie-stack))
    ;; otherwise, pop first element from trie-stack and repopulate it
    (let ((first (pop (trie--stack-store trie-stack))))
      (when first
	(trie--stack-repopulate trie-stack)
	first))))


(defun trie-stack-push (element trie-stack)
  "Push ELEMENT onto TRIE-STACK."
  (push element (trie--stack-pushed trie-stack)))


(defun trie-stack-first (trie-stack)
  "Return the first element from TRIE-STACK, without removing it
from the stack. Returns nil if the stack is empty."
  ;; if elements have been pushed onto the stack, return first of those
  (if (trie--stack-pushed trie-stack)
      (car (trie--stack-pushed trie-stack))
    ;; otherwise, return first element from trie-stack
    (car (trie--stack-store trie-stack))))


(defalias 'trie-stack-p 'trie--stack-p
  "Return t if argument is a trie-stack, nil otherwise.")


(defun trie-stack-empty-p (trie-stack)
  "Return t if TRIE-STACK is empty, nil otherwise."
  (and (null (trie--stack-store trie-stack))
       (null (trie--stack-pushed trie-stack))))



;; ----------------------------------------------------------------
;;                        Inserting data

(defun trie-insert (trie key &optional data updatefun)
  "Associate DATA with KEY in TRIE.

If KEY already exists in TRIE, then DATA replaces the existing
association, unless UPDATEFUN is supplied. Note that if DATA is
*not* supplied, this means that the existing association of KEY
will be replaced by nil.

If UPDATEFUN is supplied and KEY already exists in TRIE,
UPDATEFUN is called with two arguments: DATA and the existing
association of KEY. Its return value becomes the new association
for KEY.

Returns the new association of KEY.

Note: to avoid nasty dynamic scoping bugs, UPDATEFUN must *not*
bind any variables with names commencing \"--\"."
  (if (trie--print-form trie)
      (error "Attempted to operate on trie that is in print-form")
    ;; absurb variable names are an attempt to avoid dynamic scoping bugs
    (let ((--trie-insert--updatefun updatefun)
	  --trie-insert--old-node-flag
	  (trie-depth 0)
	  (node (trie--root trie))
	  (len (length key))
	  (i -1))
      (declare (special trie-depth))
      ;; Descend trie, adding nodes for non-existent elements of KEY. The
      ;; update function passed to `trie--insertfun' ensures that existing
      ;; nodes are left intact.
      (while (< (incf i) len)
	(setq --trie-insert--old-node-flag nil)
	(setq node (funcall (trie--insertfun trie)
			    (trie--node-subtree node)
			    (let ((trie-depth (1+ trie-depth)))
			      (trie--node-create (elt key i) trie))
			    (lambda (a b)
			      (setq --trie-insert--old-node-flag t) b)))
	(incf trie-depth))
      ;; Create or update data node.
      (setq node (funcall (trie--insertfun trie)
			  (trie--node-subtree node)
			  (trie--node-create-data data)
			  ;; if using existing data node, wrap UPDATEFUN if
			  ;; any was supplied
			  (when (and --trie-insert--old-node-flag
				     --trie-insert--updatefun)
			    (lambda (new old)
			      (setf (trie--node-data old)
				    (funcall --trie-insert--updatefun
					     (trie--node-data new)
					     (trie--node-data old)))
			      old))))
      (trie--node-data node))))  ; return new data



;; ----------------------------------------------------------------
;;                        Deleting data

(defun trie-delete (trie key &optional test)
  "Delete KEY and its associated data from TRIE.

If KEY was deleted, a cons cell containing KEY and its
association is returned. Returns nil if KEY does not exist in
TRIE.

If TEST is supplied, it should be a function that accepts two
arguments: the key being deleted, and its associated data. The
key will then only be deleted if TEST returns non-nil.

Note: to avoid nasty dynamic scoping bugs, TEST must *not* bind
any variables with names commencing \"--\"."
  (if (trie--print-form trie)
      (error "Attempted to operate on trie that is in print-form")
    (let (--trie-deleted--node
	  (--trie-delete--key key))
      (declare (special --trie-deleted--node)
	       (special --trie-delete--key))
      (trie--do-delete (trie--root trie) key test
		       (trie--deletefun trie)
		       (trie--emptyfun trie)
		       (trie--cmpfun trie))
      (when --trie-deleted--node
	(cons key (trie--node-data --trie-deleted--node))))))


(defun trie--do-delete (node --trie--do-delete--seq
			     --trie--do-delete--test
			     --trie--do-delete--deletefun
			     --trie--do-delete--emptyfun
			     --trie--do-delete--cmpfun)
  ;; Delete --TRIE--DO-DELETE--SEQ starting from trie node NODE, and return
  ;; non-nil if we deleted a node. If --TRIE--DO-DELETE--TEST is supplied, it
  ;; is called with two arguments, the key being deleted and the associated
  ;; data, and the deletion is only carried out if it returns non-nil.

  ;; The absurd argument names are to lessen the likelihood of dynamical
  ;; scoping bugs caused by a supplied function binding a variable with the
  ;; same name as one of the arguments, which would cause a nasty bug when the
  ;; lambda's (below) are called.
  (declare (special --trie-deleted--node)
	   (special --trie-delete--key))
  ;; if --TRIE--DO-DELETE--SEQ is empty, try to delete data node and return
  ;; non-nil if we did (return value of --TRIE--DO-DELETE--DELETEFUN is the
  ;; deleted data, which is always non-nil for a trie)
  (if (= (length --trie--do-delete--seq) 0)
      (setq --trie-deleted--node
	    (funcall --trie--do-delete--deletefun
		     (trie--node-subtree node)
		     (trie--node-create-dummy trie--terminator)
		     (when --trie--do-delete--test
		       (lambda (n)
			 (funcall --trie--do-delete--test
				  --trie-delete--key (trie--node-data n))))
		     nil --trie--do-delete--cmpfun))
    ;; otherwise, delete on down (return value of --TRIE--DO-DELETE--DELETEFUN
    ;; is the deleted data, which is always non-nil for a trie)
    (funcall --trie--do-delete--deletefun
	     (trie--node-subtree node)
	     (trie--node-create-dummy (elt --trie--do-delete--seq 0))
	     (lambda (n)
	       (and (trie--do-delete
		     n (trie--subseq --trie--do-delete--seq 1)
		     --trie--do-delete--test
		     --trie--do-delete--deletefun
		     --trie--do-delete--emptyfun
		     --trie--do-delete--cmpfun)
		    (funcall --trie--do-delete--emptyfun
			     (trie--node-subtree n))))
	     nil --trie--do-delete--cmpfun)))



;; ----------------------------------------------------------------
;;                       Retrieving data

(defun trie-member (trie key &optional nilflag)
  "Return the data associated with KEY in the TRIE,
or nil if KEY does not exist in TRIE.

Optional argument NILFLAG specifies a value to return instead of
nil if KEY does not exist in TRIE. This allows a non-existent KEY
to be distinguished from an element with a null association. (See
also `trie-member-p', which does this for you.)"
  (if (trie--print-form trie)
      (error "Attempted to operate on trie that is in print-form")
    ;; find node corresponding to key, then find data node, then return data
    (let (node)
      (or (and (setq node (trie--node-find trie key))
	       (trie--find-data node trie))
	  nilflag))))


(defun trie-member-p (trie key)
  "Return t if KEY exists in TRIE, nil otherwise."
  (if (trie--print-form trie)
      (error "Attempted to operate on trie that is in print-form")
    (let ((flag '(nil)))
      (not (eq flag (trie-member trie key flag))))))



;; ----------------------------------------------------------------
;;                         Completing

;; Implementation Note
;; -------------------
;; For completions ranked in anything other than lexical order, we use a
;; partial heap-sort to find the k=MAXNUM highest ranked completions among the
;; n possibile completions. This has worst-case time complexity O(n log k),
;; and is both simple and elegant. An optimal algorithm (e.g. partial
;; quick-sort where the irrelevant partition is discarded at each step) would
;; have complexity O(n + k log k), but is probably not worth the extra coding
;; effort, and would have worse space complexity unless coded to work
;; "in-place" which would be highly non-trivial. (I haven't done any
;; benchmarking, though, so feel free to do so and let me know the results!)

(defun trie-complete (trie prefix &optional rankfun maxnum reverse filter)
  "Return an alist containing all completions of PREFIX in TRIE
along with their associated data, in the order defined by
RANKFUN, defaulting to \"lexical\" order (i.e. the order defined
by the trie's comparison function). If REVERSE is non-nil, the
completions are sorted in the reverse order. If no completions
are found, return nil.

PREFIX must be a sequence (vector, list or string) containing
elements of the type used to reference data in the trie. (If
PREFIX is a string, it must be possible to apply `string' to
individual elements of the sequences stored in the trie.) The
completions returned in the alist will be sequences of the same
type as KEY. If PREFIX is a list of sequences, completions of all
sequences in the list are included in the returned alist. All
sequences in the list must be of the same type.

The optional integer argument MAXNUM limits the results to the
first MAXNUM completions. Otherwise, all completions are
returned.

If specified, RANKFUN must accept two arguments, both cons
cells. The car contains a sequence from the trie (of the same
type as PREFIX), the cdr contains its associated data. It should
return non-nil if first argument is ranked strictly higher than
the second, nil otherwise.

The FILTER argument sets a filter function for the
completions. If supplied, it is called for each possible
completion with two arguments: the completion, and its associated
data. If the filter function returns nil, the completion is not
included in the results, and does not count towards MAXNUM."

  (if (trie--print-form trie)
      (error "Attempted to operate on trie that is in print-form")

    (let (node
	  (trie--complete-accumulate
	   (if rankfun
	       (heap-create  ; heap order is inverse of rank order
		(if reverse
		    `(lambda (a b) (,rankfun a b))
		  `(lambda (a b) (not (,rankfun a b))))
		(when maxnum (1+ maxnum)))
	     (make-vector 1 nil)))
	  accumulator)

      ;; wrap prefix in a list if necessary
      ;; FIXME: the test for a list of prefixes, below, will fail if the
      ;;        PREFIX sequence is a list, and the elements of PREFIX are
      ;;        themselves lists (there might be no easy way to fully fix
      ;;        this...)
      (if (or (atom prefix) (and (listp prefix) (not (sequencep (car prefix)))))
	  (setq prefix (list prefix))
	;; sort list of prefixes if sorting completions lexically
	(when (null rankfun)
	  (setq prefix
		(sort prefix (trie-construct-sortfun
			      (trie--comparison-function trie))))))

      ;; construct function to accumulate completions
      (if rankfun
	  (setq accumulator
		(trie--complete-construct-ranked-accumulator maxnum filter))
	(setq accumulator (trie--complete-construct-accumulator maxnum filter)))

      ;; accumulate completions
      (catch 'trie-complete--done
	(mapc (lambda (pfx)
		(setq node (trie--node-find trie pfx))
		(when node
		  (trie--mapc accumulator (trie--mapfun trie) node pfx
			      (if maxnum reverse (not reverse)))))
	      prefix))

      ;; return list of completions
      (cond
       ;; extract completions from heap for ranked query
       (rankfun
	(let (completions)
	  (while (not (heap-empty trie--complete-accumulate))
	    (push (heap-delete-root trie--complete-accumulate) completions))
	  completions))
       ;; reverse result list if MAXNUM supplied
       (maxnum (nreverse (aref trie--complete-accumulate 0)))
       ;; otherwise, just return list
       (t (aref trie--complete-accumulate 0)))
      )))



(provide 'trie)

;;; trie.el ends here
