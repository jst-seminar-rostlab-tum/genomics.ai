import React from 'react';

import { Stack, Link, Tooltip } from '@mui/material';
import { Link as RouterLink } from 'react-router-dom';

// Common link with label component used in the search cards.
// eslint-disable-next-line arrow-body-style
function LabeledLink({
  label, to, content, tooltip,
}) {
  return (
    <Stack>
      {!tooltip && <div>{label}</div>}
      {tooltip && <Tooltip title={tooltip} placement="right"><div>{label}</div></Tooltip>}
      <Link
        component={RouterLink}
        to={to}
        underline="hover"
      >
        {content}
      </Link>
    </Stack>
  );
}

export default LabeledLink;
