import React from "react";
import { Tabs, Tab } from "@mui/material";

const SearchTabs = (props) => {
  return (
    <Tabs value={props.value} onChange={props.onChange}>
      <Tab label="Teams" value="teams" />
      <Tab label="Institutions" value="institutions" />
      <Tab label="Users" value="users" />
      <Tab label="Projects" value="projects" />
    </Tabs>
  );
};

export default SearchTabs;
