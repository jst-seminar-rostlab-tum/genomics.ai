import React from 'react';
import Card from '@mui/material/Card';
import Avatar from '@mui/material/Avatar';
import styles from './listCard.module.css';
import stringColor from 'shared/utils/stringColor';

function ListCard({
  imageComponent, imageURL, enforceImage, title, description, nextToTitle, trailing, onClick,
}) {
  let imageComponentUsed = imageComponent;
  if (imageURL || enforceImage) {
    if (imageComponent) {
      throw new Error('You can only have either imageComponent or imageURL.');
    }
    imageComponentUsed = (
      <Avatar
        src={imageURL}
        alt={title}
        sx={{ backgroundColor: stringColor(title || ''), width: 42, height: 42 }}
      >
        {(title || '?')[0]}
      </Avatar>
    );
  }
  return (
    <Card sx={{ border: '1px solid #c8c8c8', boxShadow: '0 0 4px rgba(111, 111, 111, 0.25)' }} elevation={0}>
      <div
        className={`${styles.cardContent}
        ${onClick ? styles.clickable : ''}`}
        onClick={onClick}
        onKeyPress={onClick}
        role="button"
        tabIndex={0}
      >
        <div className={styles.start}>
          <div className={styles.cardImageWrapper}>
            {imageComponentUsed}
          </div>
          <div className={styles.text}>
            <div className={styles.titleRow}>
              <h3 className={styles.title}>{title}</h3>
              {nextToTitle}
            </div>
            <div className={styles.description}>
              {description && (
                <p style={{
                  textOverflow: 'ellipsis', overflow: 'hidden', whiteSpace: 'nowrap', width: '60vw',
                }}
                >
                  {description}
                </p>
              )}
            </div>
          </div>
        </div>
        {trailing && (
          <div className={styles.trailing}>
            {trailing}
          </div>
        )}
      </div>
    </Card>
  );
}

export default ListCard;
