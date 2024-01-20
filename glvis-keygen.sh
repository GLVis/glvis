#!/bin/bash

# Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-443271.
#
# This file is part of the GLVis visualization tool and library. For more
# information and source code availability see https://glvis.org.
#
# GLVis is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

GLVIS_KEYS_DIR="${HOME}"/.config/glvis
NAME_COMMENT_SERVER="GLVis server key"
NAME_COMMENT_CLIENT="GLVis client key"

function have_command()
{
   command -v $1 > /dev/null 2>&1
}

function check_gpg2()
{
   if ! have_command ${GPG2}; then
      printf "\nThe required command \"${GPG2}\" was not found. Stop.\n\n"
      exit 1
   fi
}

function check_certtool()
{
   if ! have_command ${CERTTOOL}; then
      printf "\nThe required command \"${CERTTOOL}\" was not found. Stop.\n\n"
      exit 1
   fi
   if { ${CERTTOOL} --help | head -1 | grep -q GnuTLS; }; then
      # newer versions of GnuTLS certtool
      return 0
   elif { ${CERTTOOL} --version | head -1 | grep -q GnuTLS; }; then
      # older versions of GnuTLS certtool
      return 0
   else
      printf "\nThe required command \"${CERTTOOL}\" is not GnuTLS certtool."
      printf " Stop.\n\n"
      exit 1
   fi
}

function gen_keys_gpg2()
{
   local IN=""
   IN="${IN}Key-Type: RSA\n"
   IN="${IN}Key-Length: 2048\n"
   IN="${IN}Subkey-Type: RSA\n"
   IN="${IN}Subkey-Length: 2048\n"
   IN="${IN}Name-Real: "'${NAME_REAL}'"\n"
   IN="${IN}Name-Comment: "'${NAME_COMMENT}'"\n"
   IN="${IN}Name-Email: "'${NAME_EMAIL}'"\n"
   IN="${IN}Expire-Date: 0\n"
   IN="${IN}%%no-protection\n" # for gpg2 version >= 2.1
   IN="${IN}%%commit"
   IN=$(printf "${IN}")
   eval "IN=\"${IN}\""
   # printf "%s\n" "${IN}"
   local GPG2_ARGS=("--homedir" "${KEY_DIR}" "--batch" "--gen-key" "-")
   printf "%s\n" "${IN}" | ${GPG2} "${GPG2_ARGS[@]}"
}

function gen_keys_certtool()
{
   local CT="${CERTTOOL}" ROLE="$1" IN=""
   cd "${KEY_DIR}" || return 1
   # Generate user CA key
   $CT --generate-privkey --outfile ca-key.pem || return 1
   # Generate self-signed user CA certificate
   IN="cn = \"${NAME_REAL}'s GLVis ${ROLE} CA Certificate\"\n"
   IN="${IN}email = \"${NAME_EMAIL}\"\n"
   IN="${IN}expiration_days = 3651\n"
   IN="${IN}ca\n"
   IN="${IN}cert_signing_key\n"
   printf "$IN" > ca-cert.cfg
   $CT --generate-self-signed --load-privkey ca-key.pem \
       --outfile ca-cert.pem --template ca-cert.cfg \
       > ca-cert.txt 2>&1 || return 1
   rm -f ca-cert.cfg ca-cert.txt

   # Generate user key
   $CT --generate-privkey --outfile key.pem || return 1
   # Generate user certificate signed with the CA key & certificate
   IN="cn = \"${NAME_REAL}'s GLVis ${ROLE} Certificate\"\n"
   IN="${IN}email = \"${NAME_EMAIL}\"\n"
   IN="${IN}expiration_days = 3650\n"
   IN="${IN}signing_key\n"
   IN="${IN}encryption_key\n"
   if [ "${ROLE}" == "Client" ]; then
      IN="${IN}tls_www_client\n"
   else
      IN="${IN}tls_www_server\n"
   fi
   printf "$IN" > cert.cfg
   $CT --generate-certificate --load-privkey key.pem \
       --outfile cert.pem --load-ca-certificate ca-cert.pem \
       --load-ca-privkey ca-key.pem --template cert.cfg \
       > cert.txt 2>&1 || return 1
   rm -f cert.cfg cert.txt
}

function check_keys_gpg2()
{
   if [ -s "${GLVIS_KEYS_DIR}"/server/pubring.gpg ] && \
      [ -s "${GLVIS_KEYS_DIR}"/server/secring.gpg ]; then
      GEN_SERVER_KEY="NO"
   fi
   if [ -s "${GLVIS_KEYS_DIR}"/client/pubring.gpg ] && \
      [ -s "${GLVIS_KEYS_DIR}"/client/secring.gpg ]; then
      GEN_CLIENT_KEY="NO"
   fi
   if [ -s "${GLVIS_KEYS_DIR}"/client/trusted-servers.gpg ]; then
      ADD_SERVER_KEY="NO"
   fi
   if [ -s "${GLVIS_KEYS_DIR}"/server/trusted-clients.gpg ]; then
      ADD_CLIENT_KEY="NO"
   fi
}

function check_keys_certtool()
{
   if [ -s "${GLVIS_KEYS_DIR}"/server/cert.pem ] && \
      [ -s "${GLVIS_KEYS_DIR}"/server/key.pem ]; then
      GEN_SERVER_KEY="NO"
   fi
   if [ -s "${GLVIS_KEYS_DIR}"/client/cert.pem ] && \
      [ -s "${GLVIS_KEYS_DIR}"/client/key.pem ]; then
      GEN_CLIENT_KEY="NO"
   fi
   if [ -s "${GLVIS_KEYS_DIR}"/client/trusted-servers.pem ]; then
      ADD_SERVER_KEY="NO"
   fi
   if [ -s "${GLVIS_KEYS_DIR}"/server/trusted-clients.pem ]; then
      ADD_CLIENT_KEY="NO"
   fi
}

# $1 - scr dir, $2 - dest dir
function add_keys_gpg2()
{
   cd "${GLVIS_KEYS_DIR}" && \
   ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/$1 --export --armor \
           --output "$1.asc" && \
   ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/$2 --no-default-keyring \
           --keyring "trusted-${1}s.gpg" --import "$1.asc"
}

# $1 - scr dir, $2 - dest dir
add_keys_certtool()
{
   cd "${GLVIS_KEYS_DIR}" && \
   cp -fp $1/ca-cert.pem $2/trusted-${1}s.pem
}

function list_keys_gpg2()
{
   # ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/server --list-secret-keys
   if [ -s "${GLVIS_KEYS_DIR}"/server/trusted-clients.gpg ]; then
      ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/server \
              --keyring "trusted-clients.gpg" --list-public-keys
   else
      printf "Server/trusted-client keys not found.\n"
   fi
   # ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/client --list-secret-keys
   if [ -s "${GLVIS_KEYS_DIR}"/client/trusted-servers.gpg ]; then
      ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/client \
              --keyring "trusted-servers.gpg" --list-public-keys
   else
      printf "Client/trusted-server keys not found.\n"
   fi
}

function list_keys_certtool()
{
   local CT="${CERTTOOL}"
   local sc="${GLVIS_KEYS_DIR}"/server/cert.pem sn="Server certificate"
   local cc="${GLVIS_KEYS_DIR}"/client/cert.pem cn="Client certificate"
   if [ -s "$sc" ]; then
      echo "$sn:"
      echo "-------------------"
      $CT --certificate-info --infile "$sc"
   else
      printf "$sn not found.\n"
   fi
   echo
   if [ -s "$cc" ]; then
      echo "$cn:"
      echo "-------------------"
      $CT --certificate-info --infile "$cc"
   else
      printf "$cn not found.\n"
   fi
}

function read_name_email()
{
   NAME_REAL="$1"
   NAME_EMAIL="$2"
   if [ "${NAME_REAL}" == "" ]; then
      if [ -t 0 ]; then
         if [ -s "${HOME}"/.gitconfig  ] && have_command git; then
            FULL_NAME=$(git config user.name)
         elif [ $(uname -s) == "Darwin" ]; then
            FULL_NAME=$(id -P ${USER} | awk -F: '{print $8}')
         elif have_command getent; then
            FULL_NAME=$(getent passwd ${USER} | awk -F: '{print $5}')
         else
            FULL_NAME=${USER}
         fi
         read -p "Enter your name [${FULL_NAME}]: " NAME_REAL
         if [ "${NAME_REAL}" == "" ]; then
            NAME_REAL="${FULL_NAME}"
         fi
      else
         return 1
      fi
   fi
   if [ "${NAME_EMAIL}" == "" ]; then
      if [ -t 0 ]; then
         if [ -s "${HOME}"/.gitconfig  ] && have_command git; then
            FULL_EMAIL=$(git config user.email)
         else
            FULL_EMAIL="${USER}@${HOST}"
         fi
         read -p "Enter your email [${FULL_EMAIL}]: " NAME_EMAIL
         if [ "${NAME_EMAIL}" == "" ]; then
            NAME_EMAIL="${FULL_EMAIL}"
         fi
      else
         return 1
      fi
   fi
}

function print_usage()
{
   printf "Usage:\n"
   printf "   $0 {-h|--help}\n"
   printf "   $0 [<var>=<value>]... {-l|--list}\n"
   printf "   $0 [<var>=<value>]... [\"Your Name\"] [\"Your Email\"]\n"
   printf "Valid variables:\n"
   printf "    keytype={x509|gpg}      (current value: ${KEYTYPE})\n"
   printf "   certtool=<certtool-prog> (current value: ${CERTTOOL}"
   printf ", keytype: x509)\n"
   printf "       gpg2=<gpg2-prog>     (current value: ${GPG2}"
   printf ", keytype: gpg)\n"
}

function select_params()
{
   # Key generation programs: gpg2, certtool
   GPG2=${gpg2:-gpg2}
   case "$OSTYPE" in
      darwin*)
         CERTTOOL=${certtool:-gnutls-certtool}
         ;;
      *)
         CERTTOOL=${certtool:-certtool}
         ;;
   esac
   # Key type to generate: gpg or x509
   KEYTYPE=${keytype:-x509}
   case "$KEYTYPE" in
      gpg)
         keytype_prog=gpg2
         ;;
      x509)
         keytype_prog=certtool
         ;;
      *)
         printf "\nInvalid keytype: '${KEYTYPE}'. Stop.\n\n"
         exit 1
         ;;
   esac
}

while [ $# -gt 0 ]; do
   case "$1" in
      -h|--help)
         select_params
         print_usage
         exit 0
         ;;
      -l|--list)
         select_params
         check_$keytype_prog
         echo
         list_keys_$keytype_prog
         exit 0
         ;;
      -*)
         printf "Unknown option: '$1'. Stop.\n\n"
         exit 1
         ;;
      *=*)
         eval $1
         shift
         ;;
      *)
         break
         ;;
   esac
done

select_params

GEN_SERVER_KEY="YES"
GEN_CLIENT_KEY="YES"
ADD_SERVER_KEY="YES"
ADD_CLIENT_KEY="YES"
check_keys_$keytype_prog

if [ "${GEN_SERVER_KEY}" == "YES" ] || [ "${GEN_CLIENT_KEY}" == "YES" ]; then
   check_$keytype_prog
   if ! read_name_email "$@"; then
      print_usage
      exit 1
   fi
elif [ "${ADD_SERVER_KEY}" == "YES" ] || [ "${ADD_CLIENT_KEY}" == "YES" ]; then
   check_$keytype_prog
fi

echo
if [ "${GEN_SERVER_KEY}" == "NO" ]; then
   printf "Server key exists.\n\n"
else
   printf "\r---------------------------\n"
   printf "   Generating server key\n"
   printf "\r---------------------------\n"
   NAME_COMMENT="${NAME_COMMENT_SERVER}"
   KEY_DIR="${GLVIS_KEYS_DIR}"/server
   echo mkdir -p "${KEY_DIR}"
   mkdir -p "${KEY_DIR}"
   if ! gen_keys_$keytype_prog "Server"; then
      printf "\nGeneration failed. Stop.\n\n"
      exit 1
   fi
   printf "\r---------------------------\n\n"
fi

if [ "${GEN_CLIENT_KEY}" == "NO" ]; then
   printf "Client key exists.\n\n"
else
   printf "\r---------------------------\n"
   printf "   Generating client key\n"
   printf "\r---------------------------\n"
   NAME_COMMENT="${NAME_COMMENT_CLIENT}"
   KEY_DIR="${GLVIS_KEYS_DIR}"/client
   echo mkdir -p "${KEY_DIR}"
   mkdir -p "${KEY_DIR}"
   if ! gen_keys_$keytype_prog "Client"; then
      printf "\nGeneration failed. Stop.\n\n"
      exit 1
   fi
   printf "\r---------------------------\n\n"
fi

if [ "${ADD_SERVER_KEY}" == "NO" ]; then
   printf "Client's trusted-servers keyring exists.\n\n"
else
   printf "\r-------------------------------------------------------\n"
   printf " Adding the server key to the client's trusted-servers\n"
   printf "\r-------------------------------------------------------\n"
   if ! add_keys_$keytype_prog "server" "client"; then
      printf "\nOperation failed. Stop.\n\n"
      exit 1
   fi
   printf "\r-------------------------------------------------------\n\n"
fi

if [ "${ADD_CLIENT_KEY}" == "NO" ]; then
   printf "Server's trusted-clients keyring exists.\n\n"
else
   printf "\r-------------------------------------------------------\n"
   printf " Adding the client key to the server's trusted-clients\n"
   printf "\r-------------------------------------------------------\n"
   if ! add_keys_$keytype_prog "client" "server"; then
      printf "\nOperation failed. Stop.\n\n"
      exit 1
   fi
   printf "\r-------------------------------------------------------\n\n"
fi
