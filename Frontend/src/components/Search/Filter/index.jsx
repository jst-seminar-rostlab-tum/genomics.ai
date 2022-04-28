import { FormGroup, Stack } from "@mui/material";
import React from "react";
import { Route } from "react-router-dom";
import { setTypeInUrl } from "shared/utils/common/utils";
import GeneralFilter from "./GeneralFilter";
import TeamsFilter from "./TeamsFilter";

const Filter = ({ path, searchParams, updateQueryParams }) => {
  return (
    <Stack>
      <FormGroup>
        <GeneralFilter
          sortBy={searchParams.get("sortBy")}
          onChange={(param, value) => updateQueryParams(param, value)}
        />
        <Route path={setTypeInUrl(path, "teams")}>
          <TeamsFilter
            visibility={searchParams.get("visibility")}
            onChange={(param, value) => updateQueryParams(param, value)}
          />
        </Route>
      </FormGroup>
    </Stack>
  );
};

export default Filter;
