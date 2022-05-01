import React from "react";

import { Stack, Link } from "@mui/material";
import { Link as RouterLink } from "react-router-dom";

// Common link with label component used in the search cards.
const LabeledLink = (props) => {
  return (
    <Stack>
      <div>{props.label}</div>
      <Link
        component={RouterLink}
        to={props.to}
        underline="hover"
      >
        {props.content}
      </Link>
    </Stack>
  );
};

export default LabeledLink;
