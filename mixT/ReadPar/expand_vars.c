#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>

#include "readpar.h"

#ifdef HAVEREADLINE
#include <stdio.h>
#include <readline/readline.h>
#include <readline/history.h>
#endif

int RP_EXPAND_VARS = EXPANDVARS;
static int RP_EXPAND_VARS_RUN=1;

int rp_modify_par_value (char *string, char mod);

int RP_expand_variables (char *string)
{
  /* Split filename into components and then expand shell variables and then
     combine again */
  char buff[8192], value[8192];
  char *word[1024];
  char delim[8192];
  int isvar[1024]; /* 1 if word started with the $, 0 otherwise */
  int iscurlvar[1024]; /* 1 if word started with the ${, 0 otherwise */
  char c;
  char *s, *pmod;

  int i, nwords,L, j, iname, LN;
  int expandit, iscurl;
  char name[1024], nametmp[1024], prompt[256]; /* to hold variable name */

  strcpy(buff,string);
  L=strlen(buff);

  nwords = 0; j = 0; iname = 0; string[0]='\0';
  for (i=0;i<L; ) {
    if (buff[i] != '$') { /* not a variable */
      string[j]=buff[i];
      j++;
      i++;
    } else { /* this starts a variable */
      i++;
      iname = 0;
      expandit = 1;
      if (buff[i]=='{') { /* this is a {variable} */
	iscurl = 1;
	i++;
	while (i<L && buff[i] != '}') {
	  name[iname]=buff[i];
	  iname ++; i++;
	}
	if (buff[i] != '}') {
	  fprintf (stderr,"Warning: variable name not terminated in %s\n",buff);
	  expandit = 0;
	} else { /* skip '}' */
	  i++;
	}
      } else { /* this is a $variable */
	iscurl = 0;
	while (i<L && (isalnum(buff[i]) || buff[i]=='_') ) {
	  name[iname]=buff[i];
	  iname++; i++;
	}
      }
      name[iname]='\0';
      string[j]='\0';
      if (expandit) {
	/* analyze name modifiers */
	strcpy (nametmp,name);
	LN = strlen(nametmp);
	pmod = NULL;
	if (LN>2) {
	  pmod = NULL;
	  if (nametmp[LN-2]==':') {
	    switch (nametmp[LN-1]) {
	    case 'h': case 't': case 'r': case 'e':
	      pmod = nametmp+LN-1;
	      nametmp[LN-2]='\0';
	      break;
	    default:
	      pmod = NULL;
	    }
	  }
	}
	if (strncmp(nametmp,"ARGV[",5)==0 && nametmp[strlen(nametmp)-1]==']')
	  { /* ARGV[x]: get the command line arg */
	    rp_expand_vars_get_arg (nametmp,value);
	    s = value;
	  } 
	else if (strncmp(nametmp,"ASKUSER[",8)==0 && nametmp[strlen(nametmp)-1]==']')
	  {
	    strcpy (prompt,nametmp+8);
	    prompt[strlen(prompt)-1]='\0';
	    strcat (prompt,": ");
#ifdef HAVEREADLINE
	    s = readline (prompt);
	    strcpy (value,s);
#else
	    printf ("%s",prompt);
	    scanf ("%s",value);
#endif
	    s = value;
	  }
	else {
	  get_command_line_par (nametmp,value);
	  if (!defined(value)) {
	    s = getenv(nametmp);
	  } else {
	    s = value;
	  }
	}
	if (s!=NULL) {
	  if (pmod != NULL) { /* we have to modify the parameter value */
	    if (s != value) strcpy (value,s);
	    rp_modify_par_value (value,*pmod);
	    s = value;
	  }
	  strcat (string,s);
	} else { /* no par and no env under this name */
	  strcat (string,"$");
	  if (iscurl) strcat(string,"{");
	  strcat (string,name);
	  if (iscurl) strcat(string,"}");
	}
      } else {
	strcat (string,"$");
	if (iscurl) strcat(string,"{");
	strcat (string,name);
	/* this is possible only if { is not closed, so don't put the closing
	  bracket */
      }
      j = strlen(string);
    }
  }
  string[j]='\0';
}

/*************************************************************************/
/* functions to set/unset variable expansion *****************************/

int RP_q_expand_vars() {
  char *s, v[1024];
  int flag;
  if (RP_EXPAND_VARS_RUN) { /* first call not preceeded by RP_set_expand_vars*/
    RP_EXPAND_VARS_RUN = 0;
    flag = RP_EXPAND_VARS;
    if (yespar("rp_expand_vars")) {
      RP_EXPAND_VARS=1;
    } else {
      if (nopar("rp_expand_vars")) {
	RP_EXPAND_VARS=0;
      } else {
	s = getenv ("RP_EXPAND_VARS");
	if (s!=NULL) {
	  strcpy (v,s);
	  lcase (v);
	  RP_EXPAND_VARS = ( 
			    strcmp(v,"y")==0   ||
			    strcmp(v,"yes")==0 ||
			    strcmp(v,"��")==0 );
	} else {
	  RP_EXPAND_VARS = flag;
	}
      }
    }
  }
  return RP_EXPAND_VARS;
}

void rp_q_expand_vars_ (int *q) {
  *q = RP_q_expand_vars();
}

int RP_set_expand_vars (int flag) {
  RP_EXPAND_VARS_RUN = 0;
  RP_EXPAND_VARS = flag;
}

void rp_set_expand_vars_ (int *flag) {
  RP_set_expand_vars (*flag);
}


/****************************************************************/
int rp_modify_par_value (char *string, char mod) {
  char buff[8192];
  char *s, *s1;
  int retval;

  strcpy (buff,string);
  retval = 0;
  switch (mod) {

  case 'h': /* Remove a trailing  pathname  component,  leaving the head. */
  case 't': /* Remove all leading pathname components,  leaving the tail. */
    s = strrchr(buff,'/');
    if (s!=NULL) {
      if (mod=='h') {
	*s = '\0'; strcpy (string,buff);
      } else {
	strcpy (string,s+1);
      }
      retval  = 1;
    }
    return retval;
    
  case 'r': /* Remove a filename extension `.xxx', leaving  the root name. */
  case 'e': /* Remove all but the extension. */
    s = strrchr(buff,'.');
    if ( s != NULL ) {
      s1 = strrchr(buff,'/');
      if (s1 == NULL || s1 < s ) {/* there are no / or / is before the . */
	if (mod == 'r') {
	  *s = '\0'; strcpy (string,buff);
	} else {
	  strcpy (string,s+1);
	}
	retval  = 1;
      }
    }
    return retval;
    
  default: return 0;
  }
}
	       
/*************************************************************************
 Parse names like ARGV[NARG-3]
*/

int rp_expand_vars_get_arg (char *spec, char *value)
{
  char *s, *p;
  char buff[1024];
  int number, index, isnumber;
  char op;
  int wasnarg;

  s = spec+5;
  /* sanity check: expression should start with NARG or a digit */
  if ( (! isdigit(*s)) && (strncmp(s,"NARG",4)!=0) ) {
    strcpy(value,spec);
    return -1;
  }
  
  strcpy (buff,s); s=buff;
  buff[strlen(buff)-1]='\0'; /* tream trailing ] */
  op='+';
  index = 0;
  wasnarg = 0;
  while ( *s ) {
    isnumber = 0;
    if ( isspace(*s) ) {
      s++;
    } else if (strncmp(s,"NARG",4)==0) {
      s+=4;
      number = iargc()-1;
      isnumber = 1;
      wasnarg = 1;
    } else if ( *s == '+' || *s == '-' || *s == '*' || *s == '/' ) {
      if (op != '_') { /* error: last char was an operation */
	strcpy(value,spec); return -1;
      }
      op = *s;
      s++;
    } else {
      number = (int)strtol (s,&p,10);
      if ( p==s ) { /* Error: no conversion was performed */
	strcpy(value,spec); return -1;
      }
      s = p;
      isnumber = 1;
    }
    if (isnumber) {
      switch (op) {
      case '+': index+=number; break;
      case '-': index-=number; break;
      case '*': index*=number; break;
      case '/': index/=number; break;
      default:  /* not valid operation */
	strcpy(value,spec); return -1;
      }
      op = '_';
    }
  }
  if (wasnarg && index <= 0 ) {
    strcpy (value,RP_UNDEFINED);
  } else {
    getarg (index,value);
    if (value[0]=='\0') {
      strcpy (value,RP_UNDEFINED);
    }
  }
}
	
      

  
