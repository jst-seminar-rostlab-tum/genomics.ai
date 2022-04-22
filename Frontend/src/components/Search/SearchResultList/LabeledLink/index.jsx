import React from "react";

import { Stack, Link } from "@mui/material";
import { Link as RouterLink } from "react-router-dom";

const LabeledLink = (props) => {
  return (
    <Stack>
      <div>{props.label}</div>
      <Link
        component={RouterLink}
        to="documentation" // TODO integrate as prop
        underline="hover"
      >
        {props.content}
      </Link>
    </Stack>
  );
};

export default LabeledLink;
