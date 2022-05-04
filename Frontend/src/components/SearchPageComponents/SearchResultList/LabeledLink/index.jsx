import React from 'react';

import { Stack, Link } from '@mui/material';
import { Link as RouterLink } from 'react-router-dom';

// Common link with label component used in the search cards.
// eslint-disable-next-line arrow-body-style
function LabeledLink({ label, to, content }) {
  return (
    <Stack>
      <div>{label}</div>
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
