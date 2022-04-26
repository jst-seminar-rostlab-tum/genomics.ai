import TwitterIcon from '@mui/icons-material/Twitter';
import EmailIcon from '@mui/icons-material/Email';
import ContentCopyIcon from '@mui/icons-material/ContentCopy';
import { IconButton } from '@mui/material';
import React from 'react';
import { colors } from 'shared/theme/colors';

const mailIconColor = colors.neutral['600'];
const copyIconColor = colors.neutral['600'];
const twitterIconColor = '#00aced';

function ShareMenu({ projectName, url }) {
  return (
    <>
      <IconButton href={`https://twitter.com/share?ref_src=twsrc%5Etfw&url=${url}`} target="_blank">
        <TwitterIcon sx={{ color: twitterIconColor }} />
      </IconButton>
      <IconButton href={`mailto:?subject=${encodeURIComponent(projectName)}&body=${url}`}>
        <EmailIcon sx={{ color: mailIconColor }} />
      </IconButton>
      <IconButton onClick={() => navigator.clipboard.writeText(url)}>
        <ContentCopyIcon sx={{ color: copyIconColor }} />
      </IconButton>
      <script async src="https://platform.twitter.com/widgets.js" charset="utf-8" />
    </>
  );
}

export default ShareMenu;
