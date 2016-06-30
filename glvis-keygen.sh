#!/bin/bash

# Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-443271. All Rights
# reserved. See file COPYRIGHT for details.
#
# This file is part of the GLVis visualization tool and library. For more
# information and source code availability see http://glvis.org.
#
# GLVis is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.

GPG2=gpg2

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
   printf "usage:\n"
   printf "   $0 [-h|--help]\n"
   printf "   $0 [-l|--list]\n"
   printf "   $0 [\"Your Name\"] [\"Your Email\"]\n"
}

case "$1" in
   -h|--help)
      print_usage
      exit 0
      ;;
   -l|--list)
      check_gpg2
      echo
      # ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/server --list-secret-keys
      ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/server \
         --keyring "trusted-clients.gpg" --list-public-keys
      # ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/client --list-secret-keys
      ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/client \
         --keyring "trusted-servers.gpg" --list-public-keys
      exit 0
      ;;
   *)
      ;;
esac

GEN_SERVER_KEY="YES"
if [ -s "${GLVIS_KEYS_DIR}"/server/pubring.gpg ] && \
   [ -s "${GLVIS_KEYS_DIR}"/server/secring.gpg ]; then
   GEN_SERVER_KEY="NO"
fi
GEN_CLIENT_KEY="YES"
if [ -s "${GLVIS_KEYS_DIR}"/client/pubring.gpg ] && \
   [ -s "${GLVIS_KEYS_DIR}"/client/secring.gpg ]; then
   GEN_CLIENT_KEY="NO"
fi
ADD_SERVER_KEY="YES"
if [ -s "${GLVIS_KEYS_DIR}"/client/trusted-servers.gpg ]; then
   ADD_SERVER_KEY="NO"
fi
ADD_CLIENT_KEY="YES"
if [ -s "${GLVIS_KEYS_DIR}"/server/trusted-clients.gpg ]; then
   ADD_CLIENT_KEY="NO"
fi

if [ "${GEN_SERVER_KEY}" == "YES" ] || [ "${GEN_CLIENT_KEY}" == "YES" ]; then
   check_gpg2
   if ! read_name_email "$@"; then
      print_usage
      exit 1
   fi
elif [ "${ADD_SERVER_KEY}" == "YES" ] || [ "${ADD_CLIENT_KEY}" == "YES" ]; then
   check_gpg2
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
   if ! gen_keys_gpg2; then
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
   if ! gen_keys_gpg2; then
      printf "\nGeneration failed. Stop.\n\n"
      exit 1
   fi
   printf "\r---------------------------\n\n"
fi

cd "${GLVIS_KEYS_DIR}"

if [ "${ADD_SERVER_KEY}" == "NO" ]; then
   printf "Client's trusted-servers keyring exists.\n\n"
else
   printf "\r-------------------------------------------------------\n"
   printf " Adding the server key to the client's trusted-servers\n"
   printf "\r-------------------------------------------------------\n"
   ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/server --export --armor\
      --output "server.asc" && \
   ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/client --no-default-keyring\
      --keyring "trusted-servers.gpg" --import "server.asc"
   if [ $? -ne 0 ]; then
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
   ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/client --export --armor\
      --output "client.asc" && \
   ${GPG2} --homedir "${GLVIS_KEYS_DIR}"/server --no-default-keyring\
      --keyring "trusted-clients.gpg" --import "client.asc"
   if [ $? -ne 0 ]; then
      printf "\nOperation failed. Stop.\n\n"
      exit 1
   fi
   printf "\r-------------------------------------------------------\n\n"
fi
