import TwitterIcon from '@mui/icons-material/Twitter';
import EmailIcon from '@mui/icons-material/Email';
import { IconButton } from '@mui/material';
import React from 'react';
import { colors } from 'shared/theme/colors';

const mailIconColor = colors.neutral['600'];
const twitterIconColor = '#00aced';

function ShareMenu({ projectName, url }) {
  return (
    <>
      <IconButton aria-label="mail" href={`mailto:?subject=${encodeURIComponent(projectName)}&body=${url}`}>
        <EmailIcon sx={{ color: mailIconColor }} />
      </IconButton>
      <IconButton href={`https://twitter.com/share?ref_src=twsrc%5Etfw&url=${url}`} target="_blank">
        <TwitterIcon sx={{ color: twitterIconColor }} />
      </IconButton>
      <script async src="https://platform.twitter.com/widgets.js" charset="utf-8" />
    </>
  );
}

export default ShareMenu;
