#!/bin/bash

dotnet tool restore
dotnet paket restore
dotnet build src/Uncodium.Eigensystems.sln --configuration Release
